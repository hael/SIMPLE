# solve3D chunk chains: explicit stream emulation and automatic balanced chunks

Implementation note, 2026-10-10; extended for `solve3D_chunks` on 2026-10-11.
Both workflows are implemented but uncompiled and not yet run; implementation
step 5 (validation, section 9) remains open, so this single living note stays
under `planned`. The legacy `emulate_solve3D_stream` decisions were reviewed
by Hans on 2026-10-10 (section 10). Written against master `81b4e6741` plus
the uncommitted implementation described here. Companion policies:
`doc/policies/3D/solve3D_policy.md`, `doc/policies/3D/solve3D_addon_policy.md`.

## 1. Purpose

Stream task 7 (`simple_stream_stage_solve3D`) runs one `solve3D` on the first
particles the pool publishes and then one `solve3D_addon` per growth of the
pool. Its run time and the quality of its maps depend on a handful of settings
(base size, cohort size, `nstates`, `nstages`, `lpstop`, `nsample`,
`rec_backend`, `nparts`/`nthr`) that cannot be tuned in a live session: a
session cannot be replayed, and its timing is entangled with the upstream
tasks.

There are now two offline `simple_exec` entry points over the same
`solve3D`/`solve3D_addon` chain.

`emulate_solve3D_stream` is the explicit experimental control. It preserves
its original interface and behavior: the caller fixes the base and add-on
sizes, and can choose whether regressed add-ons are rolled back.

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

`solve3D_chunks` is the production-facing automatic policy. The caller gives
the project and the stream-compatible `solve3D` inputs, but no base size,
add-on size, `nsample`, or rollback policy:

```text
simple_exec prg=solve3D_chunks projfile=<project> <stream-compatible solve3D inputs ...>
```

It samples a fixed total of `NSAMPLE_TOTAL_CHUNK=5000` particles per
iteration whatever `nstates` is, derives the run count from that total,
deals the active particles to the runs class by class, runs partition 0 as the base, and runs every later partition once
as an add-on. The sampling, the split and its refusal are the contract in
section 3.2.

## 2. What the stream does, and what the emulation reproduces

This comparison remains the contract of `emulate_solve3D_stream`; adding
`solve3D_chunks` does not change it.

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

`solve3D_chunks` reuses the same prepared-source, step-project, nested
commander, and reporting mechanics. It is not a closer imitation of the
stream cadence: its purpose is a deterministic automatic analysis of all
active particles in balanced, capped searched cohorts. In particular, it
does not use the stream's row order, whole-stack rounding, cadence gate, or
rollback.

## 3. Partitioning

Both workflows count the rows whose input `ptcl2D` state is > 0. The
emulation takes them in row order (the stream's order: pool stacks in import
order); `solve3D_chunks` deals them by 2D class (3.2). Inactive rows are
assigned to no partition and stay inactive in every step. Let `Nactive` be
the number of active rows.

### 3.1 Explicit legacy emulation

The `emulate_solve3D_stream` partition rule is unchanged:

```text
base       = active particles 1 .. nptcls_base
rem        = Nactive - nptcls_base
nchunks    = max(1, rem / nptcls_addon)          (integer division)
chunk k    = nptcls_addon particles, k = 1 .. nchunks-1
chunk last = rem - (nchunks-1)*nptcls_addon      (in [nptcls_addon, 2*nptcls_addon),
                                                   or rem when rem < nptcls_addon)
```

Example: `Nactive=10000, nptcls_base=3000, nptcls_addon=2000` gives
3000 | 2000 | 2000 | 3000. When the remainder is below `nptcls_addon`, there
is one add-on on that remainder; the base is never enlarged because its size
is the experiment's control variable.

Before launching a child run, the emulation refuses
`nptcls_base < 5*nstates`, `nptcls_addon < 5*nstates`,
`rem < 5*nstates`, and `nptcls_base >= Nactive`. Consequently the legacy
emulation still requires at least one add-on and has no automatic 30,000
particle cap.

### 3.2 Automatic balanced chunks

`solve3D_chunks` has no caller-selected base or add-on count. Its constants
and derived policy are:

```text
NSAMPLE_TOTAL_CHUNK = 5000                          particles sampled per iteration, all states together
nsample             = NSAMPLE_TOTAL_CHUNK / nstates (rounded down)
NMIN                = 2 * NSAMPLE_TOTAL_CHUNK = 10000
NMAX                = 6 * NSAMPLE_TOTAL_CHUNK = 30000
nruns               = ceil(Nactive / NMAX)
```

solve3D samples `nsample * nstates` particles per iteration
(`update_frac = nsample * nstates / N`). That policy is not changed. The
workflow divides the total by `nstates` instead, so `nsample * nstates` is
5000 rounded down to a multiple of `nstates` (4998 for three states), and the
number of states influences neither the number of particles sampled nor the
run sizes. The update fraction of every run stays between 1/6 and 1/2. Each
state gets about `5000 / nstates` freshly searched particles per iteration.

The `Nactive` particles are dealt to exactly `nruns` runs, class by class:
ordered by their `ptcl2D` class label (row order within a class, unlabelled
particles first), the `i`-th goes to run `mod(i-1, nruns)`. Run sizes differ
by at most one: if `Nactive = q*nruns + r`, the first `r` runs contain `q+1`
active particles and the rest contain `q`. Every class is split over the
runs to within one particle, so every run, the base included, has the same
2D class mix, and a cohort is a fair sample of the base's population. The
runs are not contiguous in the rows. Without a usable 2D solution (`ptcl2D`
virgin or no `cls2D`) every label counts as 0 and the deal follows row
order. Run 0 is the base; it is not assigned a separate `nptcls_base`. Runs
1 through `nruns-1` are the add-on cohorts. The full contract and its
reasons are in `doc/policies/3D/solve3D_chunks_policy.md`.

`nruns=1` is valid, so a data set of 10,000 to 30,000 active particles runs
the base only. Before any child `solve3D` run the workflow refuses
`nstates < 1` and `Nactive < NMIN`. Nothing else can fail: with `NMAX = 3 *
NMIN`, every part of a split into two or more runs holds at least
`NMAX / 2 = 15000` particles.

The 30,000-particle maximum applies to the newly searched cohort of each
`solve3D` or `solve3D_addon` invocation. It is not a cap on the frozen union:
after each adopted add-on, the frozen project contains the cumulative union
and can grow past 30,000 until the final union contains all `Nactive`
particles.

The explicit and balanced partitioners both return one partition label per
project row and are unit-testable without running either reconstruction
commander.

## 4. Steps

### 4.1 Directory layout

The workflow's `params%new` (mkdir=yes) creates the execution directory `E`
and copies the input project into it (`E/<project>.simple`); the workflow
prepares `source.simple` from that copy. All steps run from `E`:

```text
E/source.simple                     input copy, prepared once (4.2), never run on
E/base.simple                       base step project (rows of partitions >= 1 deselected)
E/1_solve3D/base.simple             base result, registers its run manifest
E/addon_01.simple                   add-on 1 project; replaced by the add-on's published union
E/2_solve3D_addon/                  add-on 1 run directory (frozen copy, report, maps)
E/addon_02.simple, E/3_solve3D_addon/ ...
E/emulate_solve3D_stream_report.txt legacy emulation step table (section 6)
E/solve3D_chunks_report.txt          automatic-chunk step table (section 6)
```

One invocation writes only the report corresponding to its entry point. Both
reports are rewritten after every completed step. `solve3D_chunks` also
publishes its final result over `E/<project>.simple`, where `simple_exec`
records the job (the emulation does not).

### 4.2 Preparing the source copy

For either workflow, once on `source.simple`: drop everything 3D and every
run registration
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

Every step in either workflow executes the existing commanders (`commander_solve3D`,
`commander_solve3D_addon`) commander-style, within the code environment, as
the library's workflows do (`exec_solve3D_cavgs_conditional_restarts` runs
`solve3D_cavgs` this way). Each workflow's executor records `E` before a
step and returns to it after: `solve3D` leaves the working directory in its
run directory, and `solve3D_addon` publishes by absolute path without
returning.

### 4.4 Verdicts, rollback, and fixed-cohort adoption

After add-on k, both workflows read `solve3D_addon_report.txt` from its run
directory and retain the per-state verdicts in their workflow report.

The legacy `emulate_solve3D_stream` adoption rule is unchanged:

- `rollback=yes` (default, the stream's rule): a REGRESSED state leaves the
  frozen project unchanged. Partition k stays selected in add-on k+1's
  project, so its cohort is partitions k and k+1, the closest deterministic
  analogue of the stream's retry with a grown cohort. A regressed last
  add-on leaves the previous result as the final one.
- `rollback=no`: every add-on result is adopted.

Either way, a regressed emulation run's directory and published project stay
for inspection, and the report flags the step.

`solve3D_chunks` has no rollback input and always adopts every add-on,
including one whose diagnostic verdict is REGRESSED. This preserves the
precomputed fixed cohorts: each balanced partition is newly searched exactly
once, the frozen union advances by exactly that partition, and no rejected
partition is merged into a later add-on. A regressed add-on is logged with
the states it affects and stays visible in `solve3D_chunks_report.txt`; the
verdict does not change the plan.

## 5. Command line

The UI entries share the private `simple_ui_solve3D` routine
`add_solve3D_base_inputs`. `solve3D` asks for the full list;
`emulate_solve3D_stream` and `solve3D_chunks` ask for the stream-compatible
list, which omits `vol1`, `cavg_ini`, `cavg_ini_ext`, and `state`. A separate
logical argument controls whether the public `nsample` input is present:
it remains present for `solve3D` and the legacy emulation, and is absent for
`solve3D_chunks`.

### 5.1 Explicit legacy emulation

| Key | Status | Routed to |
| --- | --- | --- |
| `nptcls_base` (existing key) | required | partitioning |
| `nptcls_addon` (new key) | required | partitioning |
| `rollback` (new key, yes\|no, default yes) | optional | the emulation |
| `addon_diag` (existing key) | optional | every add-on |
| `solve3D` inputs except those below | as in `solve3D` | the base |
| compute and add-on keys the `solve3D_addon` UI accepts (`nparts`, `nthr`, `nsample`, `overlap`, `pcg_solvent_check`) | as given | the base and every add-on |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state` | refused | none: input volumes, class-average starts and state continuation are not the stream's base, and their class averages would not match the base subset |

This interface and its defaults are preserved. In particular, `nsample` is
still caller-controlled, `rollback=yes` is still the default, and explicit
base and add-on sizes are still required.

### 5.2 Automatic balanced chunks

| Key | Status | Routed to |
| --- | --- | --- |
| `projfile` | required | prepared source and every step project |
| `addon_diag` | optional, default no | every add-on |
| stream-compatible `solve3D` inputs | as in the emulation | the base; manifest-recorded policy reaches add-ons |
| `nsample` | not accepted from the caller | set internally to `NSAMPLE_TOTAL_CHUNK / nstates` for the base and replayed from its manifest by add-ons |
| `nptcls_base`, `nptcls_addon`, `rollback` | not accepted | derived balanced plan and unconditional adoption replace them |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state`, `projfile_frozen` | not accepted | none |

`solve3D_chunks` checks its registered UI allowlist, so hidden manual chunk
controls are refused rather than silently ignored. The only chunk-policy
knobs are the internal constants in section 3.2.

### 5.3 Shared routing and defaults

The add-on forwarding rule is mechanical, so it follows the UIs. An entry
command-line key goes to an add-on when the `solve3D_addon` UI accepts it and
the base manifest does not record it (`manifest_records_input`). Recorded
keys (`balance`, `nclust`, `mskdiam`, `maxits_pcg`, `maxits_ml`, ...) reach
the add-ons through the base run's manifest. Workflow-only keys stay behind.
`nsample` is not a recorded input, but `replay` sets it, with `nstates`, from
the base run's sampling record, so the add-ons sample the base run's total.
`solve3D_chunks` therefore leaves its internal `nsample` off the add-on
command line.

Defaults are `solve3D`'s, `balance` included: `cavg` when the project
carries class averages (kept in the source copy, 4.2), `none` when it has no
usable 2D solution (`set_balance_default` in `exec_solve3D`; the `class`
branch needs `vol1`, which is refused here). Neither workflow replaces this
default. An explicit `balance` is honoured.

## 6. Report

The shared report type writes `emulate_solve3D_stream_report.txt` for the
legacy entry point and `solve3D_chunks_report.txt` for the automatic entry
point. Each has one row per completed step and is rewritten after every step,
so a run that stops part-way keeps what it measured.

| Column | Source |
| --- | --- |
| step, kind (base \| addon), chunk index | the owning workflow executor and partition plan |
| frozen, added (newly searched cohort), union size | the adopted chain and partition sizes |
| last stage run, stage `lp` and `box_crop` | the step's manifest (`get_last_stage`, `get_stage`) |
| wall time (s), and wall time per 1000 cohort particles | `system_clock` around the commander call |
| adopted \| rolled back | section 4.4; automatic chunks are always adopted |
| per state: FSC=0.143 and FSC=0.5 resolution of the union | the result's registered FSC (`get_fsc`, `get_resolution`) |
| per state, add-ons: verdict, FSC shell shift, map correlation, cohort-only res (with `addon_diag=yes`) | `solve3D_addon_report` getters (`get_verdict`, `get_dshell`, `get_corr`, `get_cohort_res0143`) |

The legacy header retains its original command line, source project, selected
count, requested base and add-on sizes, number of add-on steps, `nstates`, and
rollback policy. The `solve3D_chunks` header instead records the selected
count, number of runs, `nstates`, the derived `nsample` with the total it
samples per iteration (`nsample * nstates`), and every planned run size.

## 7. Architecture

| Piece | Location | Content |
| --- | --- | --- |
| partition and run-plan policy | `src/main/solve/simple_solve3D_stream_emulation.f90` | legacy `partition_selected`; automatic `chunk_nsample`, `derive_chunk_run_count` and `partition_selected_balanced`; `NSAMPLE_TOTAL_CHUNK`, `NPTCLS_CHUNK_MIN` and `NPTCLS_CHUNK_MAX` |
| shared step-project utilities | the same solve module | `prepare_emulation_source` and `select_emulation_step`, used by both workflows |
| shared report utilities | the same solve module | `emulation_step`, the backward-compatible `emulation_report%new`, the chunks-specific `new_chunks`, and the two report filenames |
| partition/report tester | `src/main/solve/simple_solve3D_stream_emulation_tester.f90` | both partition policies, balanced feasibility boundaries, inactive rows, source/step projects, and shared report serialization |
| commanders | renamed `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90` | public `commander_emulate_solve3D_stream` and `commander_solve3D_chunks`; each owns its validation, partition choice, child-command sequence, key exclusions, and adoption policy, while sharing only leaf helpers for step-project writing, accepted argument copying, and statistics |
| UI | `src/main/ui/simple/simple_ui_solve3D.f90` | both entries and shared `add_solve3D_base_inputs`, including the flag that hides `nsample` from `solve3D_chunks` |
| dispatch | `src/main/exec/simple_exec_solve3D.f90` | both `emulate_solve3D_stream` and `solve3D_chunks` cases, importing the renamed commander module |
| keys | `simple_parameters.f90`, `simple_parameters_parse.f90` | legacy `nptcls_addon` and `rollback`; `solve3D_chunks` adds no public parameter keys |
| policy | `doc/policies/3D/solve3D_addon_policy.md` | a short section naming the emulation and its differences from the stream (section 2) |

The commander module owns orchestration but no reconstruction algorithm. Its
two independent executors call `commander_solve3D` and
`commander_solve3D_addon` and read their published projects, manifests, and
add-on reports. Only policy-neutral leaf routines are shared between those
executors. The solve module owns the small deterministic run-plan,
project-selection, and text-report utilities shared by the two workflows.
The underlying solve3D commanders and the live stream stage are unchanged;
replacing the stream's `cap_rows` with either offline partitioner would be a
separate change.

## 8. Risks and checks before trusting the numbers

- **Row growth is not exercised.** With every row present from the start the
  superset relation never sees appended rows (`check_appended_stacks`), and
  the layout digest never changes row count. Neither affects timing or maps.
  The stream's append path stays covered by the snapshot generator test.
- **Partitioning versus whole stacks.** The stream's caps round down to whole
  stacks. Both offline workflows partition active particles exactly; the
  automatic workflow additionally balances run sizes and every class's share
  to within one. A stack-aware mode would be a different policy because it
  cannot in general preserve that balance.
- **Dealt cohorts read every stack.** A `solve3D_chunks` run reads its
  particles from all stacks rather than from a block of stacks: the same
  amount of data from more files.
- **The cap is a search-work cap, not a memory cap for the union.** Every
  automatic cohort is at most 30,000 particles, but each add-on also reads
  and reconstructs with the growing frozen union. Late-step cost can
  therefore grow even though newly searched work remains capped.
- **Fewer sampled particles per state.** The total sampled per iteration is
  fixed, so with `nstates` states each state gets about `5000 / nstates`
  freshly searched particles per iteration. Multi-state runs may need more
  iterations to turn over every state's particles.
- **Automatic chunks favor fixed experiments over stream rollback.** A
  REGRESSED add-on is reported but adopted. Users who need the stream's
  retry behavior retain `emulate_solve3D_stream rollback=yes`.
- **Load balance.** In the emulation, deselected rows of later partitions sit
  at the end of the row range in every step, as in a capped stream run. If
  partitioning by row range leaves late parts nearly empty, the timing
  measures that imbalance, which the stream has too. `solve3D_chunks` spreads
  every step's active rows over the whole row range.
- **Look-ahead in the 2D solution** (4.2).

## 9. Implementation steps

1. Preserve the legacy keys, explicit partitioner, source preparation,
   step-project construction, rollback path, and report coverage.
2. Add the automatic constants, the `nsample` derivation, the run count,
   the active-balanced partitioner, the report run-size vector, and boundary
   tests.
3. Rename the commander module to `simple_commanders_solve3D_chunks`, add
   the second public commander, preserve the legacy executor as an independent
   workflow, and share only small policy-neutral helpers with the new executor.
4. Register and dispatch `solve3D_chunks`; share the stream-compatible UI
   inputs while withholding manual chunk controls and public `nsample`.
5. Open validation: run the registered unit suite and exercise both entry
   points on representative projects, including a legacy rollback case, an
   automatic base-only case, a split with a one-particle balance remainder,
   the 30,000 boundary, multi-state runs (the logged and replayed
   `nsample * nstates` total), report contents, and final union membership.
   Compilation and runtime validation have not yet been performed.
6. After validation, record the evidence here and move this note out of
   `planned`; keep this file as the single implementation record for both
   workflows.

## 10. Decisions

Decisions 1-9 are the legacy `emulate_solve3D_stream` decisions reviewed on
2026-10-10 and remain unchanged:

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

The `solve3D_chunks` extension adds these decisions:

10. Chunk analysis derives all sizing from a fixed sampling total,
    `NSAMPLE_TOTAL_CHUNK=5000` particles per iteration, with
    `NMIN=2*NSAMPLE_TOTAL_CHUNK=10000` and `NMAX=6*NSAMPLE_TOTAL_CHUNK=30000`;
    none is a command-line input. The number of states influences neither the
    total sampled nor the run sizes: the workflow sets
    `nsample=NSAMPLE_TOTAL_CHUNK/nstates`, and solve3D's own policy
    (`nsample * nstates` sampled per iteration) stays unchanged (Hans,
    2026-10-11).
11. Use the fewest capped runs, `nruns=ceil(Nactive/NMAX)`, and deal the
    active particles to them class by class, so run sizes and every class's
    share differ by at most one. Run 0 is the base in the same deal, not a
    separately sized special case. Without 2D classes the deal follows row
    order (Hans, 2026-10-11).
12. Permit a base-only plan and refuse fewer than `NMIN` active particles;
    no other plan can fail.
13. Always adopt automatic add-ons, even when their quality verdict is
    REGRESSED, so every fixed balanced cohort is searched exactly once.
    Rollback remains available only in the explicit legacy emulation.
14. Interpret 30,000 as the maximum newly searched cohort in one child run,
    not as a maximum on the cumulative frozen union.

## 11. As implemented

- The renamed module
  `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90` exports
  both `commander_emulate_solve3D_stream` and
  `commander_solve3D_chunks`. Each type-bound `execute` routine contains its
  own complete orchestration; there is no mode switch or shared top-level
  driver.
- The legacy entry remains developer-visible and still requires
  `nptcls_base` and `nptcls_addon`, accepts caller `nsample`, and defaults to
  `rollback=yes`. Its explicit partition and rejected-cohort merge behavior
  are unchanged.
- The advanced `solve3D_chunks` entry accepts the same stream-compatible
  solve3D settings but no `nptcls_base`, `nptcls_addon`, `rollback`, or
  `nsample`. It checks the UI allowlist, reads `nstates` from the command
  line, injects `nsample=chunk_nsample(nstates)` before `params%new`, passes
  the `ptcl2D` class labels (zeros without a usable 2D solution) to the
  balanced partitioner, permits zero add-ons, and sets rollback false so
  every result advances the frozen union. It logs and reports `nsample` with
  the total `nsample * nstates`. It refuses `nstages` below
  `solve3D_symsrch_stage()` before anything is written, logs every regressed
  add-on, and at the end publishes the final result over its execution
  directory's project copy (`publish_solve3D_result`, the pattern of
  `solve3D_addon`'s publication: temporary file renamed over the target,
  manifest registered by absolute path).
- `src/main/solve/simple_solve3D_stream_emulation.f90` remains the shared
  utility module. It now contains the legacy and balanced partitioners,
  `chunk_nsample`, `derive_chunk_run_count`, the fixed constants, common
  source preparation and step selection, and the `emulation_step` and
  `emulation_report` types.
  The original report initializer and legacy output remain intact; a separate
  chunks initializer records the balanced run-size vector and how the runs
  were dealt (by 2D class, with the class count, or in row order) in
  `solve3D_chunks_report.txt`. The balanced partitioner deals by a counting
  sort on the class label; it no longer uses `split_nobjs_active`.
- The input project is copied into the workflow execution directory by its
  single `params%new`; the prepared source is `E/source.simple`. Every step
  project is rewritten from it with absolute `projfile` entries. Both child
  commanders run in-process, and each workflow returns to `E` after each one
  with `simple_chdir`.
- A step's run directory is found from its registered manifest. The add-on
  publishes the union over `E/addon_NN.simple` and registers that manifest by
  absolute path; the shared collector reads the add-on report and FSC data
  there. Per-state FSC=0.143 and FSC=0.5 resolutions use the registered FSC,
  with `refine3D_fsc_fname` in the run directory as fallback.
- Each workflow explicitly owns its forwarding exclusions. A small shared
  leaf helper copies a key only when `solve3D_addon` accepts it and the base
  manifest does not record it. Legacy workflow keys are removed from the
  legacy base; `addon_diag` is removed from both bases; automatic `nsample`
  stays out of add-on overrides because it is replayed from the manifest.
- `simple_solve3D_stream_emulation_tester.f90` covers both partition
  policies, including the automatic constants, the `nsample` derivation for
  1 to 20 states, the cap-derived run count and its bounds over a sweep of
  particle counts, the deal in row order without labels, the deal by class
  on a two-session fixture (equal sizes, every class split to within one,
  every class in the base), base-only bounds, and the pre-existing
  source/project/report behavior. These tests are written
  but have not yet been compiled or run.
