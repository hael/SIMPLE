# emulate_solve3D_stream: the stream's solve3D/solve3D_addon cycle, offline

Implementation note, 2026-10-10. Planned; nothing here is implemented.
Decisions reviewed by Hans on 2026-10-10 (section 10).
Written against master `d205b3189` plus the uncommitted add-on stage-policy
change (registration pass + the base run's last stage). Companion policies:
`doc/policies/3D/solve3D_policy.md`, `doc/policies/3D/solve3D_addon_policy.md`.

## 1. Purpose

Stream task 7 (`simple_stream_stage_solve3D`) runs one `solve3D` on the first
particles the pool publishes and then one `solve3D_addon` per growth of the
pool. Its run time and the quality of its maps depend on a handful of settings
(base size, cohort size, `nstates`, `nstages`, `lpstop`, `nsample`,
`rec_backend`, `nparts`/`nthr`) that cannot be tuned in a live session: a
session cannot be replayed, and its timing is entangled with the upstream
tasks.

`emulate_solve3D_stream` replays the cycle on an existing project, under
`simple_exec`, with the cohort schedule fixed by two numbers:

```text
simple_exec prg=emulate_solve3D_stream projfile=<project> nptcls_base=<n> nptcls_addon=<m> \
            <solve3D inputs ...>
```

The first `nptcls_base` selected particles are the base run; the rest are added
in add-on runs of `nptcls_addon` particles each, the last partition merged
with the previous when it falls below `nptcls_addon`. One invocation emulates
one setting; a parameter sweep is several invocations. The result is the same
chain of projects the stream would publish plus one report that puts the
steps' timing and quality side by side.

## 2. What the stream does, and what the emulation reproduces

From `simple_stream_stage_solve3D.f90` (start_solve3D, start_addon,
cap_rows, roll_back_addon) and `simple_commanders_stream_p07_solve3D_multistate.f90`:

| Stream behaviour | Emulation |
| --- | --- |
| One growing row set; particles not in a run are deselected (state 0) in the job's project only (`set_queued_states`) | same: every step's project is the full row set with the rows of later partitions deselected |
| Base: `solve3D` with `balance=none`, `pgrp=c1`, `nstates`, `nstages`, `mskdiam`, `nparts3D`, `nthr3D` | base: `solve3D` with the emulation's command line (section 5); `balance` follows `solve3D`'s default rule (`cavg` with class averages, `none` without) |
| Add-on: `solve3D_addon projfile=<grown> projfile_frozen=<latest result>`, compute keys only; everything else replayed from the manifest | same, through the same commander |
| Cohort capped at `nptcls_addon_max`, whole stacks in order (`cap_rows`); cadence gate `max(5*nstates, 0.10*nfrozen)` (`next_job`) | fixed partitions of `nptcls_addon` particles (section 3); no cadence gate, no whole-stack rounding |
| A REGRESSED verdict rolls the run back: the previous result stays the frozen base and the cohort is retried with more particles (`roll_back_addon`, `retry_cohort`) | `rollback=yes` (default): the previous result stays frozen and the rejected partition joins the next (section 4.4) |
| Row growth from new pool stacks (`n > nf` in the superset relation) | not reproduced: every step holds all rows (`n == nf`); see section 8 |

The emulation does not reproduce the upstream tasks (the pool's 2D, its
selection model, the sieve's `cavg_ini_ext` start) or the final multistate
`refine3D` the stream runs on `pool_final`.

## 3. Partitioning

The particles are the rows whose `ptcl2D` state is > 0 in the input project,
in row order (the stream's order: pool stacks in import order). Let `N` be
their count.

```text
base      = active particles 1 .. nptcls_base
rem       = N - nptcls_base
nchunks   = max(1, rem / nptcls_addon)          (integer division)
chunk k   = nptcls_addon particles, k = 1 .. nchunks-1
chunk last= rem - (nchunks-1)*nptcls_addon      (in [nptcls_addon, 2*nptcls_addon), or rem when rem < nptcls_addon)
```

Example: `N=10000, nptcls_base=3000, nptcls_addon=2000` gives 3000 | 2000 |
2000 | 3000.

A partition owns the rows from its first active particle up to the row before
the next partition's first active particle; inactive rows inside that range
stay inactive in every step. When the remainder is below `nptcls_addon` there
is a single add-on on the remainder: the base is never enlarged, since its
size is the experiment's control variable.

Refusals, before anything is written: `nptcls_base < 5*nstates` (the
stream's `MIN_PTCLS_PER_STATE`), `nptcls_addon < 5*nstates`
(`MIN_COHORT_STATE_POP` of `simple_project_superset`), `rem < 5*nstates`, and
`nptcls_base >= N`.

The partitioning is one pure function (input: the active flags in row order
and the three integers; output: a partition label per row, 0 for inactive
rows), so it is unit-testable without a project.

## 4. Steps

### 4.1 Directory layout

`simple_exec` (mkdir=yes) creates the execution directory `E` and copies the
input project into it as `source.simple`. All steps run from `E`:

```text
E/source.simple                     input copy, prepared once (4.2), never run on
E/base.simple                       base step project (rows of partitions >= 1 deselected)
E/1_solve3D/base.simple             base result, registers its run manifest
E/addon_01.simple                   add-on 1 project; replaced by the add-on's published union
E/2_solve3D_addon/                  add-on 1 run directory (frozen copy, report, maps)
E/addon_02.simple, E/3_solve3D_addon/ ...
E/emulate_solve3D_stream_report.txt the step table (section 6), rewritten after every step
```

### 4.2 Preparing the source copy

Once, on `source.simple`: drop everything 3D and every run registration
(`cls3D`, the 3D entries of `out`, the `ptcl3D` alignment, the `projinfo`
entries `sigma2_state`, `solve3D_manifest`, `solve3D_run_id`), as the snapshot
generator does. The 2D solution (`cls2D`, the class-average and `frc2D`
entries of `out`) is kept when present, so `balance=cavg` and the
FRC-planned ladder remain available. They are the full data set's classes,
a mild look-ahead the stream does not have (its first publication's FRCs come
from the particles published so far). Without a usable 2D solution the base
runs under `balance=none` with the `mskdiam` ladder, as `solve3D` now does
(`solve3D_policy.md`, "Without a 2D solution").

### 4.3 Base and add-on k

Every step project is a copy of `source.simple` with its `ptcl2D` and
`ptcl3D` states set to the source's state for the rows of partitions
0..k and to 0 for the rest. It is rebuilt from the source each time: the
frozen particles' poses and labels come from the frozen project, never from
the current one (`project_superset%restore`), so the step projects need
nothing but the selection.

- **Base.** `solve3D projfile=E/base.simple mkdir=yes` plus the base command
  line (section 5). Its result `E/1_solve3D/base.simple` is the first frozen
  project. The base must write a run manifest: a single-state run stopped by
  `nstages` at stage 3 or later now does (final reconstruction), stages 1-2
  are refused up front.
- **Add-on k.** `solve3D_addon projfile=E/addon_kk.simple
  projfile_frozen=<current frozen> mkdir=yes` plus the forwarded keys. The
  commander publishes the union over `E/addon_kk.simple`, which then registers
  the add-on's manifest and is the frozen project of add-on k+1. Each add-on
  runs the registration pass and the base run's last stage.

Every step executes the existing commanders (`commander_solve3D`,
`commander_solve3D_addon`) commander-style, within the code environment, as
the library's workflows do (`exec_solve3D_cavgs_conditional_restarts` runs
`solve3D_cavgs` this way). The emulation records `E` before a step and
returns to it after: `solve3D` leaves the working directory in its run
directory, and `solve3D_addon` publishes by absolute path without
returning.

### 4.4 Verdicts and rollback

After add-on k the emulation reads `solve3D_addon_report.txt` from its run
directory (`solve3D_addon_report%read`, as `read_addon_verdict` does).

- `rollback=yes` (default, the stream's rule): a REGRESSED state leaves the
  frozen project unchanged. Partition k stays selected in add-on k+1's
  project, so its cohort is partitions k and k+1, the closest deterministic
  analogue of the stream's retry with a grown cohort. A regressed last
  add-on leaves the previous result as the final one.
- `rollback=no`: every add-on result is adopted.

Either way the regressed run's directory and published project stay for
inspection, and the report flags the step.

## 5. Command line

The UI entry carries the `solve3D` inputs and the emulation's own. To keep
the two in step, the `solve3D` input list moves into a private subroutine of
`simple_ui_solve3D` (`add_solve3D_base_inputs`) that `new_solve3D` and
`new_emulate_solve3D_stream` both call.

| Key | Status | Routed to |
| --- | --- | --- |
| `nptcls_base` (existing key) | required | partitioning |
| `nptcls_addon` (new key) | required | partitioning |
| `rollback` (new key, yes\|no, default yes) | optional | the emulation |
| `addon_diag` (existing key) | optional | every add-on |
| `solve3D` inputs except those below | as in `solve3D` | the base |
| compute and add-on keys the `solve3D_addon` UI accepts (`nparts`, `nthr`, `nsample`, `overlap`, `pcg_solvent_check`) | as given | the base and every add-on |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state` | refused | none: input volumes, class-average starts and state continuation are not the stream's base, and their class averages would not match the base subset |

The forwarding rule is mechanical, so it follows the UIs. A key on the
emulation's command line goes to an add-on when the `solve3D_addon` UI
accepts it and the manifest does not record it (`manifest_records_input`).
The recorded keys (`balance`, `nclust`, `mskdiam`, `maxits_pcg`, `maxits_ml`,
...) reach the add-ons through the base run's manifest, so forwarding them
would only log overrides equal to the replayed value. The emulation's own
keys are stripped from the base command line.

Defaults are `solve3D`'s, `balance` included: `cavg` when the project
carries class averages (kept in the source copy, 4.2), `none` when it has no
usable 2D solution (`set_balance_default` in `exec_solve3D`; the `class`
branch needs `vol1`, which is refused here). The emulation sets no default of
its own. An explicit `balance` is honoured.

## 6. Report

`emulate_solve3D_stream_report.txt`, one row per step, rewritten after every
step so a run that stops part-way keeps what it measured:

| Column | Source |
| --- | --- |
| step, kind (base \| addon), run directory, result project | the emulation |
| partition range (first/last active particle), `ncohort`, `nfrozen`, union size | the partitioning; for add-ons the counts the add-on logs (`project_superset`) |
| last stage run, stage `lp` and `box_crop` | the step's manifest (`get_last_stage`, `get_stage`) |
| wall time (s), and wall time per 1000 cohort particles | `system_clock` around the commander call |
| per state: FSC=0.143 and FSC=0.5 resolution of the union | the result's registered FSC (`get_fsc`, `get_resolution`) |
| per state, add-ons: verdict, FSC shell shift, map correlation, cohort-only res (with `addon_diag=yes`) | `solve3D_addon_report` getters (`get_verdict`, `get_dshell`, `get_corr`, `get_cohort_res0143`) |
| adopted \| rolled back | section 4.4 |

The header records the command line, the partition sizes and the source
project, so sweep results can be collated by script.

## 7. Architecture

| Piece | Location | Content |
| --- | --- | --- |
| partitioning | new `src/main/solve/simple_solve3D_stream_emulation.f90` | the pure partition function (section 3) and the report type; no commander or I/O dependencies beyond the report file |
| partition tester | new `src/main/solve/simple_solve3D_stream_emulation_tester.f90` | the merge rule, a remainder below one chunk, inactive rows inside partitions, every refusal |
| step projects | the same module, one routine on `sp_project` | the source preparation (4.2) and a step project from the source and a partition bound (4.3) |
| commander | new `src/main/commanders/simple/simple_commanders_solve3D_stream_emulation.f90` | `commander_emulate_solve3D_stream`: validation, the step loop, key routing, verdicts, report; no solve3D logic |
| UI | `src/main/ui/simple/simple_ui_solve3D.f90` | `new_emulate_solve3D_stream`, and the shared `add_solve3D_base_inputs` |
| dispatch | `src/main/exec/simple_exec_solve3D.f90` | `case('emulate_solve3D_stream')` |
| keys | `simple_parameters.f90`, `simple_parameters_parse.f90` | `nptcls_addon` (int), `rollback` (char) |
| policy | `doc/policies/3D/solve3D_addon_policy.md` | a short section naming the emulation and its differences from the stream (section 2) |

Boundaries: the emulation calls the two commanders and reads their published
outputs (projects, manifests, add-on reports); it changes neither commander
and adds no branch to `exec_solve3D`. The stream stage is not touched; if
the emulation's partitioning later replaces `cap_rows`, that is a separate
change.

## 8. Risks and checks before trusting the numbers

- **Row growth is not exercised.** With every row present from the start the
  superset relation never sees appended rows (`check_appended_stacks`), and
  the layout digest never changes row count. Neither affects timing or maps.
  The stream's append path stays covered by the snapshot generator test.
- **Partitioning versus whole stacks.** The stream's caps round down to whole
  stacks. Exact partitions make the cohort size the controlled variable. A
  `chunk_mode=stacks` option can follow if the rounding turns out to matter.
- **Load balance.** Deselected rows of later partitions sit at the end of the
  row range in every step, as in a capped stream run. If partitioning by row
  range leaves late parts nearly empty, the timing measures that imbalance,
  which the stream has too.
- **Look-ahead in the 2D solution** (4.2).

## 9. Implementation steps

1. Keys `nptcls_addon`, `rollback`. The partition function and its tester.
2. Source preparation and step-project construction, with a tester on a
   synthetic project (states only).
3. `add_solve3D_base_inputs` refactor (no behaviour change; the UI JSON of
   `solve3D` must be byte-identical before and after), then the new UI entry.
4. The commander: validation and refusals, the base step, the add-on loop,
   key routing, rollback, the report. Dispatch in `simple_exec_solve3D`.
5. A validation run on the snapshot test data: base of 1500, add-ons of
   500, single-state `nstages=5` and multi-state `nstates=3`.
6. Policy section; the test inventory entry.

## 10. Decisions (reviewed 2026-10-10)

1. Masked steps on one row set, instead of physical snapshots (the snapshot
   generator's copied stacks): faithful to the stream's capped runs, no image
   I/O, only the append path is not exercised (section 8).
2. Counts are of selected particles, not rows (section 3). Agreed.
3. A remainder below one chunk becomes a single add-on, never merged into the
   base. Agreed.
4. `rollback=yes` by default, with the rejected partition joining the next.
5. `balance` follows `solve3D`'s default: `cavg` with class averages, `none`
   without. Agreed.
6. The source project's 2D classes and FRCs are kept, despite their
   look-ahead. Agreed for this experiment.
7. Add-ons get only forwarded keys (section 5); add-on settings of their own
   (for example their own `nsample`) are not in the first version. Agreed.
8. The steps execute the `solve3D` and `solve3D_addon` commanders
   commander-style within the code environment, as elsewhere in the library;
   no subprocess route. Agreed.
9. Later, optional: a batch reference `solve3D` on all particles with the
   same settings, with the final union compared against it (map correlation
   and FSC, using the add-on report's comparison routines).
