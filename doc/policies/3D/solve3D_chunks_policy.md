# solve3D_chunks Policy

This document records the policy for `solve3D_chunks`. It solves a 3D map
from all active particles of a project in balanced chunks: a `solve3D` run on
the first chunk, then one `solve3D_addon` run per later chunk, each frozen on
the result before it. The contracts of the two runs are in
[solve3D_policy.md](solve3D_policy.md) and
[solve3D_addon_policy.md](solve3D_addon_policy.md). The design record shared
with `emulate_solve3D_stream` is
[emulate_solve3D_stream.md](../../implementation_notes/planned/emulate_solve3D_stream.md).

Status (2026-10-11): implemented, not yet compiled or run.

## 1. Scope

`simple_exec prg=solve3D_chunks projfile=<project> <solve3D inputs>` takes a
project with a selection (`ptcl2D` state > 0) and returns a 3D solution of
every selected particle in the project of its execution directory. No run
searches more than 30000 particles. The caller chooses no chunk sizes. The
workflow derives them from the number of selected particles alone.

The workflow orchestrates only. It does not change the sampling, staging,
reconstruction or validation of `solve3D` and `solve3D_addon`. It runs both
unchanged, in the same process, and reads what they publish.

It is not the stream. It has no cadence gate, no whole-stack rounding and no
rollback. Those belong to the stream and to `emulate_solve3D_stream`, which
replays the stream's cycle with caller-chosen chunk sizes.

## 2. Command Line

| Key | Status |
| --- | --- |
| `projfile` | required |
| `addon_diag` | optional, default `no`; goes to every add-on |
| the `solve3D` inputs a stream base run has (`nstates`, `nstages`, `pgrp`, `mskdiam`, `balance`, `nclust`, `lpstart`, `lpstop`, `rec_backend`, ...) | as in `solve3D`; they go to the base run; `nstages` must be at least 3 (section 9) |
| compute keys (`nparts`, `nthr`) and the execution-environment keys | accepted; they go to the base run and every add-on |
| `nsample` | not accepted; derived (section 3) |
| `nptcls_base`, `nptcls_addon`, `rollback` | not accepted; the plan and the adoption rule replace them |
| `vol1`, `cavg_ini`, `cavg_ini_ext`, `state`, `projfile_frozen` | not accepted |

The workflow checks every key against its own UI entry and refuses anything
else by name. Nothing is ignored silently.

`mkdir` must be `yes` or absent. The workflow always works in its own
numbered execution directory, on a copy of the input project. The input
project is never written; the copy receives the result (section 7).

Defaults are `solve3D`'s, `balance` included: `cavg` when the project carries
class averages, `none` when it has no usable 2D solution. The workflow sets no
default of its own except `nsample`.

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

It deals the selected particles to `R` parts, class by class:

1. Order the selected particles by their 2D class label (`ptcl2D` class),
   keeping row order within a class. Particles without a label (0) form one
   group that comes first.
2. Give the `i`-th particle of that order to part `mod(i-1, R)`.

The parts are not contiguous in the rows. Within each part's project the
rows keep their order.

This gives:

- **Equal sizes.** Part sizes differ by at most one. If `N = q*R + r`, the
  first `r` parts hold `q+1` particles and the rest hold `q`.
- **Equal class mix.** Each class is split over the parts to within one
  particle, so every part, the base included, has the same 2D class
  composition. Groups of similar class averages (`balance=cavg`) are then
  balanced too.
- **Determinism.** No random numbers are drawn; the same project always gives
  the same parts.

Part 0 is the base run. Parts 1 to `R-1` are the add-on cohorts, in order.
The base gets no special size or content.

The reason is that the base is the reference of the whole chain: every
add-on is aligned against the frozen maps of the runs before it. Cut in row
(import) order, the base would be the earliest-imported particles. Session
drift, and views or classes that are rare early in the data, would carry
into the base map, and a later cohort could differ from the base as a
population, so that its `REGRESSED` verdict could not be told from a
misalignment. Dealt by class, every cohort is a fair sample of the same
population as the base. The split uses only the 2D classes, never maps or
poses.

Without a usable 2D solution (`ptcl2D` never classified, or no `cls2D`) all
labels count as 0 and the particles are dealt in row order. That still
spreads session drift over the runs, but not view imbalance.

Unselected rows belong to no part and stay unselected in every run.
`solve3D_addon` supports cohorts interleaved with the frozen rows; its
distributed partitions follow the active rows
([solve3D_addon_policy.md](solve3D_addon_policy.md), section 9). Each run reads its particles from all stacks rather than from a
block of stacks: the same amount of data from more files.

With 10000 to 30000 selected particles there is one part, and only the base
run is executed.

The only refusal is fewer than 10000 selected particles. Nothing else can
fail: with two or more runs, `N > 30000 * (R-1)`, so the smallest part holds
at least 15000 particles.

Examples:

| N | runs | part sizes |
| --- | --- | --- |
| 9999 | refused | — |
| 10000 | 1 | 10000 |
| 30000 | 1 | 30000 |
| 30001 | 2 | 15001, 15000 |
| 60001 | 3 | 20001, 20000, 20000 |
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
SAME`) and stays visible in the report. For the stream's rollback, use
`emulate_solve3D_stream rollback=yes`.

## 6. Key Routing

- **Base run:** the entry command line with `prg=solve3D`, the base step
  project, the derived `nsample`, and `addon_diag` removed.
- **Add-ons:** a fresh command line with `projfile`, `projfile_frozen` and
  `mkdir=yes`. An entry key is added when the `solve3D_addon` UI accepts it
  and the base manifest does not record it as an input (`addon_diag`,
  `nparts`, `nthr`, `overlap`, `pcg_solvent_check`). Recorded inputs
  (`balance`, `nclust`, `mskdiam`, `maxits_pcg`, `maxits_ml`, ...) reach the
  add-ons through the manifest.
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
E/2_solve3D_addon/                  add-on 1 run directory: frozen copy, maps, manifest, solve3D_addon_report.txt
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
- the planned size of every part, and whether the parts were dealt by 2D
  class (with the number of classes) or in row order;
- one row per run: kind, part, frozen and added particles, union size, last
  stage with its low-pass limit and box, wall time, and seconds per 1000
  added particles (add-ons);
- per state and run: the FSC=0.143 and FSC=0.5 resolutions of the result;
  for add-ons also the verdict, the FSC shell shift, the map correlation
  against the frozen map, and with `addon_diag=yes` the cohort-only
  resolution.

## 8. Cost

The cap bounds the searched work of each run, not the union. Every add-on
also holds the frozen particles of all earlier parts and reconstructs the
whole union at its end. Late add-ons therefore cost more than early ones,
even though each searches the same number of particles.

`solve3D_addon` warns when a cohort is below 5 % of the frozen population.
With equal parts that happens from the 21st add-on on, that is with more
than 630000 selected particles.

## 9. Refusals

Before anything is written:

- a key that is not an input of `solve3D_chunks`;
- `mkdir` other than `yes`;
- `nstates < 1`;
- `nstages` below 3, the symmetry-search stage. Such a base run would write
  no manifest (one state), or `solve3D_addon` would refuse to continue it
  (several states). In practice `nstages` is at least 5; the refusal only
  makes a mistake fail at once instead of after the base run.

After the project copy, before any step project or child run:

- the copy did not land in the workflow's own execution directory;
- no particles in `ptcl2D`, or `ptcl2D` and `ptcl3D` of different size;
- fewer than 10000 selected particles.

The child runs keep their own refusals.

## 10. Tests

- Unit sub-suite `solve3D stream emulation` in `unit_project`: the sampling
  constants; `nsample` for 1 to 20 states (the total stays within
  `nstates - 1` of 5000); the run count at its bounds and over a sweep of
  particle counts (fewest runs, every part between 10000 and 30000); the deal
  in row order without labels, with unselected rows skipped; the deal by
  class on two sessions in row order (equal run sizes, every class split to
  within one particle, every class in the base); base-only plans; the
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
| `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90` (`exec_solve3D_chunks`, `publish_solve3D_result`) | validation, `nsample`, the class labels for the deal, the run sequence, key routing, statistics, the regression log, the report file, publication |
| `src/main/solve/simple_solve3D_stream_emulation.f90` (`chunk_nsample`, `derive_chunk_run_count`, `partition_selected_balanced`, `prepare_emulation_source`, `select_emulation_step`, `emulation_report`) | the sampling constants, the run count, the deal by class, the source and step projects, the report format |
| `src/main/solve/simple_solve3D_stream_emulation_tester.f90` | the unit tests |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the child runs, unchanged |
| `src/main/solve/simple_solve3D_manifest.f90` (`replay`, `read_registered`, `register`) | `nsample`, `nstates` and the recorded inputs reaching the add-ons; the manifest registration of the published project |
| `src/main/solve/simple_solve3D_utils.f90` (`solve3D_symsrch_stage`) | the lowest supported `nstages` |
