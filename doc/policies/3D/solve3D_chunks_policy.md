# solve3D_chunks Policy

This document records the policy for `solve3D_chunks`. It solves a 3D map
from all active particles of a project in balanced chunks: a `solve3D` run on
the first chunk, then one `solve3D_addon` run per later chunk, each frozen on
the result before it. The contracts of the two runs are in
[solve3D_policy.md](solve3D_policy.md) and
[solve3D_addon_policy.md](solve3D_addon_policy.md). The design record is
[emulate_solve3D_stream_completed.md](../../implementation_notes/completed/emulate_solve3D_stream_completed.md);
the explicit-size emulator it also describes, `emulate_solve3D_stream`, was
removed on 2026-10-14.

## 1. Scope

`simple_exec prg=solve3D_chunks projfile=<project> <solve3D inputs>` takes a
project with a selection (`ptcl2D` state > 0) and returns a 3D solution of
every selected particle in the project of its execution directory. No run
searches more than 30000 particles. The caller chooses no chunk sizes. The
workflow derives them from the number of selected particles alone.

The workflow orchestrates only. It does not change the sampling, staging,
reconstruction or validation of `solve3D` and `solve3D_addon`. It runs both
unchanged, in the same process, and reads what they publish.

Its runs are consecutive whole stacks in arrival order, with sizes chosen
automatically within the limits below. It has no cadence gate and no rollback.

## 2. Command Line

| Key | Status |
| --- | --- |
| `projfile` | required |
| `addon_diag` | optional, default `yes`; goes to every add-on |
| `balance` | optional, default as in `solve3D`: `cavg` with a usable 2D solution, `none` without |
| `filt_mode` | optional, default `nonuniform_lpset` |
| the other `solve3D` inputs a stream base run has (`nstates`, `nstages`, `pgrp`, `mskdiam`, `nclust`, `lpstart`, `lpstop`, ...) | as in `solve3D`; they go to the base run; `nstages` must reach the symmetry-search stage (section 9) |
| compute keys (`nparts`, `nthr`) and the execution-environment keys | accepted; they go to the base run and every add-on |
| `nsample` | not accepted; derived (section 3) |
| `rec_backend`, `maxits_pcg`, `maxits_ml`, `pcg_solvent`, `pcg_solvent_lambda`, `pcg_solvent_check` | not accepted; reconstruction is always `rec_backend=gridding` |
| `nptcls_base` | not accepted; the plan replaces it |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state`, `projfile_frozen` | not accepted |

The workflow checks every key against its own UI entry and refuses anything
else by name. Nothing is ignored silently.

`mkdir` must be `yes` or absent. The workflow always works in its own
numbered execution directory, on a copy of the input project. The input
project is never written; the copy receives the result (section 7).

The workflow sets its own defaults before the base run, and the add-ons
inherit them through the base manifest (`addon_diag` is forwarded):

- `addon_diag=yes`: every add-on also reconstructs its cohort alone, which
  cross-validates the cohort's alignment against the frozen map.
- `filt_mode=nonuniform_lpset`: nonuniform filtering, with the NU frontier
  fed into the LP-set matching run from the NU filtering stage on. `solve3D`
  promotes its own default, `nonuniform`, to the same mode at that stage.
- `rec_backend=gridding`, fixed.

Every other default is `solve3D`'s, `balance` included: `cavg` when the
project has a usable 2D solution, `none` when it has none.

## 3. Sampling

`solve3D` samples `nsample * nstates` particles per iteration. Its update
fraction is `nsample * nstates / N`, where `N` is the run's number of searched
particles. That rule is `solve3D`'s and is not changed here.

`solve3D_chunks` fixes the total instead. Every run samples 5000 particles per
iteration, whatever the number of states. The workflow sets

`nsample = 5000 / nstates` (rounded down)

on the base run, and `solve3D` multiplies it back by `nstates`:

| nstates | nsample | sampled per iteration |
| --- | --- | --- |
| 1 | 5000 | 5000 |
| 2 | 2500 | 5000 |
| 3 | 1666 | 4998 |
| 4 | 1250 | 5000 |

Rounding down loses fewer than `nstates` particles, and the total never
exceeds 5000. So the number of states influences neither how many particles
are sampled nor the chunk sizes (section 4).

The base run records `nsample` and `nstates` in its manifest. Every add-on
replays both, so every add-on samples the same total from its own cohort.

With a fixed total, each state gets about `5000 / nstates` freshly searched
particles per iteration. Multi-state runs therefore turn over each state's
particles more slowly than single-state runs.

Every run is fractional. `solve3D` forces full sampling only when
`nsample / N > 0.9`; with `nsample <= 5000` and `N >= 10000` that ratio is at
most 0.5.

## 4. Partitioning

Let `N` be the number of selected particles. Every run searches between
10000 and 30000 new particles: two to six times the sampling total. The
update fraction of every run therefore lies between 1/6 and 1/2.

The workflow uses the fewest runs the cap allows:

`R = ceiling(N / 30000)`

It cuts the selected particles, in row order, into `R` runs of consecutive
whole stacks: the earliest-imported stacks first, no micrograph split.

1. The balanced targets: if `N = q*R + r`, the first `r` runs would hold
   `q+1` particles and the rest `q`.
2. Cut `k` goes to the stack boundary nearest its target (ties to the earlier
   boundary), among the boundaries that keep run `k-1` and the remaining runs
   within 10000 to 30000.
3. When no such boundary exists (in practice only with stacks of thousands of
   particles), the cut goes inside a stack, at the target moved into the
   limits.

Run sizes therefore differ from the balanced share by at most about one stack,
and every run holds 10000 to 30000 particles. The log and the report give the
number of cuts that fell inside a stack. The split is deterministic and uses
neither class labels nor maps.

Run 0 is the base run. Runs 1 to `R-1` are the add-on cohorts, in arrival
order.

The runs keep the order of arrival:

- **The base is the start of the session.** Session drift, and views or
  classes that are rare early in the data, carry into the base map.
- **View imbalance within a run** is left to the sampler: `balance=cavg`
  evens out the views among the particles a run has. Views a cohort lacks come
  from the frozen term of its add-on.
- **Each run reads its own block of stack files.**
  The particles of a run sit together on disk.

Unselected rows belong to no run and stay unselected in every run.

With 10000 to 30000 selected particles there is one run, and only the base
run is executed.

The only refusal is fewer than 10000 selected particles. Nothing else can
fail: with two or more runs, `N > 30000 * (R-1)`, so the balanced share is at
least 15000 and a cut within the limits always exists.

Examples, with stacks of 50 particles:

| N | runs | run sizes |
| --- | --- | --- |
| 9999 | refused | — |
| 10000 | 1 | 10000 |
| 30000 | 1 | 30000 |
| 30001 | 2 | 15000, 15001 |
| 100000 | 4 | 25000 each |

## 5. Run Sequence

Let `E` be the execution directory.

1. **Source copy.** The input project is copied into `E`. From it the
   workflow writes `E/source.simple` with everything 3D removed: the 3D
   alignment, update counts, `cls3D`, state maps and FSCs, state weights, and
   the run registrations (`sigma2_state`, `solve3D_manifest`,
   `solve3D_run_id`). The 2D solution stays. `ptcl3D` states follow the
   `ptcl2D` selection.
2. **Step projects.** The project of step `k` is the source with parts
   `k+1` and later deselected in `ptcl2D` and `ptcl3D`. Each step therefore
   holds parts 0 to `k`. All steps share one row set; rows are never added or
   removed.
3. **Base run.** `solve3D` runs on `E/base.simple` (part 0) with the entry
   command line and the derived `nsample`. Its result is
   `E/1_solve3D/base.simple`.
4. **Add-on `k`** (`k = 1 .. R-1`). `solve3D_addon` runs on
   `E/addon_kk.simple` (parts 0 to `k`), frozen on the previous result. Its
   cohort is part `k`. It publishes the union over `E/addon_kk.simple`, which
   becomes the frozen project of add-on `k+1`.
5. After every run the workflow returns to `E`, reads the run's manifest,
   FSCs and add-on report, and rewrites the report (section 7).
6. **Publication.** The final result replaces the project copy in `E`
   (section 7).

**Every add-on is adopted**, including one whose report holds a `REGRESSED`
verdict. Each part is searched exactly once, and the frozen union grows by
exactly that part. A rejected part is never merged into a later cohort, which
could break the 30000 cap. A regression is logged with the states it affects
(`SOLVE3D_CHUNKS ADDON STEP k REGRESSED IN STATE(S): ...; ADOPTED ALL THE
SAME`) and stays visible in the report.

## 6. Key Routing

- **Base run:** the entry command line with `prg=solve3D`, the base step
  project, the derived `nsample`, and `addon_diag` removed.
- **Add-ons:** a fresh command line with `projfile`, `projfile_frozen` and
  `mkdir=yes`. An entry key is added when the `solve3D_addon` UI accepts it
  and the base manifest does not record it as an input (`addon_diag`,
  `nparts`, `nthr`, `overlap`). Recorded inputs (`balance`, `nclust`,
  `mskdiam`, `filt_mode`, ...) reach the add-ons through the manifest; the
  backend does not, since every add-on reconstructs by gridding.
- **`nsample`** is never put on an add-on command line. It is not a recorded
  input, but the manifest's `replay` sets it, together with `nstates`, from
  the base run's sampling record.

## 7. Outputs

```text
E/<project>.simple                  copy of the input project, replaced by the final result at the end
E/source.simple                     the 3D-free source (section 5)
E/base.simple                       base step project
E/1_solve3D/base.simple             base result
E/addon_01.simple                   add-on 1 step project, replaced by its published union
E/2_solve3D_addon/                  add-on 1 run directory: frozen copy, maps, manifest, solve3D_addon_report.txt, BENCH/
E/addon_02.simple, E/3_solve3D_addon/, ...
E/solve3D_chunks_report.txt         the report
```

The final result is the last add-on's project, or `E/1_solve3D/base.simple`
when only the base ran. When the chain has finished, the workflow publishes
it over `E/<project>.simple`, the project in which `simple_exec` records the
job, the way `solve3D_addon` publishes its own result: it writes the project
beside the target under a temporary name and renames it over the target, and
it registers the run manifest by absolute path. The published project
carries the final union's solution and its manifest registration, as
`solve3D_addon`'s own published result does. A run that fails part-way
publishes nothing; the copy then still holds the input project. The log names the
final result and the publication (`SOLVE3D_CHUNKS FINAL RESULT`,
`SOLVE3D_CHUNKS PUBLISHED THE FINAL RESULT TO`).

`solve3D_chunks_report.txt` is rewritten after every run, so a run that stops
part-way keeps what it measured. It holds:

- the source project and the command line, with the derived `nsample`;
- the number of selected particles, runs and states;
- `nsample` and the total it samples per iteration (`nsample * nstates`);
- the planned size of every run, and the number of run cuts that fell
  inside a stack;
- one row per run: kind, part, frozen and added particles, union size, last
  stage with its low-pass limit and box, wall time, and seconds per 1000
  added particles (add-ons);
- per state and run: the FSC=0.143 and FSC=0.5 resolutions of the result;
  for add-ons also the verdict, the FSC shell shift, the map correlation
  against the frozen map, and the cohort-only resolution (`addon_diag=yes`,
  the default).

## 8. Cost

The cap bounds the searched work of each run, and the add-on passes bound its
reads: the first add-on reads the base run's particles once per box it needs
them at (the stage box and the native box) and publishes the first carried
sums; every later add-on adopts the previous run's carried sums and reads only
its own particles, in the stage iterations and in the incremental final
reconstruction (solve3D_addon policy, sections 5 and 8). A late add-on
therefore reads no more particle images than an early one with the same
cohort; what still grows with the union is the size of the carried sums it
copies and the sigma2 state it streams.

Every `solve3D` and `solve3D_addon` run writes `BENCH/SOLVE3D_BENCH.txt` in its
run directory (SIMPLE benchmarking, `L_BENCH_GLOB`; refine3D policy, section
11.1): the times of setup, the
starting state (for an add-on, the frozen accumulation), the add-on
registration pass, the refine3D stages, the final reconstruction and the
epilogue (for an add-on, its report and the cohort-only map). Read across the
chain, they separate the capped search from the union passes.

`solve3D_addon` warns when a cohort is below 5 % of the frozen population.
With equal parts that happens from the 21st add-on on, that is with more
than 630000 selected particles.

## 9. Refusals

Before anything is written:

- a key that is not an input of `solve3D_chunks`;
- `mkdir` other than `yes`;
- `nstates < 1`;
- `nstages` below the symmetry-search stage (3; 4 when `pgrp_start`
  differs from `pgrp`). Such a base run would write
  no manifest (one state), or `solve3D_addon` would refuse to continue it
  (several states). In practice `nstages` is at least 5; the refusal only
  makes a mistake fail at once instead of after the base run.

After the project copy, before any step project or child run:

- the copy did not land in the workflow's own execution directory;
- no particles in `ptcl2D`, or `ptcl2D` and `ptcl3D` of different size;
- fewer than 10000 selected particles.

The child runs keep their own refusals.

## 10. Tests

- Unit sub-suite `solve3D chunks plan` in `unit_project`: the sampling
  constants; `nsample` for 1 to 20 states (the total stays within
  `nstates - 1` of 5000); the run count at its bounds and over a sweep of
  particle counts (fewest runs, every run between 10000 and 30000); runs of
  whole stacks with unselected rows inside, sizes within one stack of the
  balanced share, a cut inside a single large stack; base-only plans; the
  refusals; the report header.
- No workflow gate exists yet. Validation still has to cover a base-only run,
  a split with a one-particle remainder, the 30000 boundary, a multi-state
  run (the logged and replayed total), the report, the regression log line,
  the early `nstages` refusal, and the published project (its solution and
  its manifest registration).

## 11. Code Map

| File | Owns |
| --- | --- |
| `src/main/ui/simple/simple_ui_solve3D.f90` (`new_solve3D_chunks`, `add_solve3D_base_inputs`) | the command-line contract |
| `src/main/exec/simple_exec_solve3D.f90` | the program entry |
| `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90` (`exec_solve3D_chunks`, `publish_solve3D_result`) | validation, `nsample`, the defaults, the stack indices for the cuts, the run sequence, key routing, statistics, the regression log, the report file, publication |
| `src/main/solve/simple_solve3D_chunks_plan.f90` (`chunk_nsample`, `derive_chunk_run_count`, `partition_selected_balanced`, `prepare_chunks_source`, `select_chunks_step`, `chunks_report`) | the sampling constants, the run count, the whole-stack cuts, the source and step projects, the report format |
| `src/main/solve/simple_solve3D_chunks_plan_tester.f90` | the unit tests (suite `solve3d_chunks_plan`) |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` (`write_solve3D_benchmark`) | the child runs; their `BENCH/SOLVE3D_BENCH.txt` |
| `src/main/solve/simple_solve3D_manifest.f90` (`replay`, `read_registered`, `register`) | `nsample`, `nstates` and the recorded inputs reaching the add-ons; the manifest registration of the published project |
| `src/main/solve/simple_solve3D_controller.f90` (`solve3D_symsrch_stage`) | the lowest supported `nstages` |
