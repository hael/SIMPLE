# solve3D_addon Policy

This document records the current policy for `solve3D_addon`, which extends
a completed 3D solution with particles it did not contain. The base contracts
are in [solve3D_policy.md](solve3D_policy.md) and
[refine3D_policy.md](refine3D_policy.md). The design history, the review
record and the decisions behind this policy are in
[abinitio3D_addon_mode_proposal.md](../../implementation_notes/completed/abinitio3D_addon_mode_proposal.md).

## 1. Scope

`solve3D_addon projfile=<current> projfile_frozen=<solution>` takes a
completed solution (the frozen project) and a current project that holds its
particles and more. The particles of the frozen solution (the frozen
particles) keep their poses and state labels and are never searched; they
contribute their accumulated signal to every reconstruction. The other active
particles of the current project (the cohort) are searched against the frozen
maps, from stage 3 of the frozen run's planned ladder to its last stage. The
result is one ordinary project with every particle posed.

Two uses are supported:

- **A more permissive selection of a fixed project.** The frozen solution was
  computed on a harsh selection; the add-on adds what a more permissive
  selection of the same particles contains. The two projects have the same
  rows and differ in their selection.
- **A growing stream.** Each update appends newly classified particles to the
  previous update's particles and is frozen on the previous update's output,
  so it searches only the new particles (section 12).

It is not a refinement of the frozen particles, not a class-average route
(the add-on aligns particles only and never runs `solve3D_cavgs`) and not
an ini3D phase.

## 2. Command Line

The add-on runs with the frozen run's settings; its command line carries only
what governs compute, convergence and diagnostics. The contract is the
program's UI entry (`new_solve3D_addon` in `simple_ui_solve3D.f90`),
checked by the wrapper through `ui_program%accepts`
([ui_layer_policy.md](../ui_layer_policy.md)); there is no second key list:

| Group | Inputs the UI entry declares |
| --- | --- |
| Projects | `projfile_frozen` (`projfile` and `mkdir` as for every project program) |
| Compute | `nparts`, `nthr` |
| Sampling and convergence | `nsample` (default: the frozen run's effective value), `overlap` (default 0.95 at stage 3) |
| PCG solve budget and checks | `maxits_pcg`, `maxits_ml`, `pcg_solvent_check` |
| Diagnostics | `addon_diag` |

The execution environment every program accepts from its launcher passes
through unchanged (`UI_ENVIRONMENT_KEYS`: the queue system, NICE, and the
stream's `worker_server` and `worker_priority`). Every accepted key is
forwarded as given.

Refusal is by key, so an inherited value cannot be confirmed silently. A key
the frozen run's manifest records as one of its inputs
(`manifest_records_input`: the solution, reconstruction and search policy, and
the entry routes as provenance) is refused as set by the frozen run; any
other key the UI entry does not accept is refused as not an input.
Ordinary `solve3D` and `solve3D_cavgs` refuse `projfile_frozen` and
`addon_diag`.

`center=no` is forced: a re-centred reference would map a shift onto cohort
particles that the frozen accumulators never receive.

## 3. The Frozen Input

The run manifest is the only route into the add-on. `exec_solve3D` writes
`solve3D_manifest.txt` at the end of every completed run (final
reconstruction or `refine3D_states` handoff) and registers it in `projinfo`;
an early-stopped run writes none and is not a frozen input. The frozen project
must be an eligible `solve3D` or `solve3D_addon` output whose manifest
validates against the project that registered it:

- the run identifier registered in `projinfo` matches the manifest's;
- row count, particle layout digest, stack table and optics/CTF parameters
  are the ones the manifest recorded;
- every state map registered in the project is the map the run produced
  (digest);
- the committed residual sigma2 state the manifest records exists unchanged;
- the current project's native box and sampling equal the frozen solution's;
- no stage box of the inherited ladder exceeds the native box.

The manifest reader keeps no backwards compatibility (release 4): a field
or input key it does not know, such as those of removed features
(`ptcl_src`, `objfun_den`, `conical_fsc`, `inpl_cont`), makes the manifest
unreadable, so a solution from an older build cannot be extended.

The frozen project is never written. The add-on works on a copy of it in the
run directory (`frozen/`, with its sigma2 state), whose sigma2 file must stay
byte-equal to the frozen run's through every frozen accumulation.

## 4. Particle Identity and Membership

The two projects share one particle index space: row `i` names the same
particle image in both wherever both hold a row `i`. They may differ in size.

- Every row both projects hold must name, in `ptcl2D` and in `ptcl3D`, the
  same stack file, image index, stack box and sampling.
- Rows the current project appends past the frozen project's last row are
  cohort candidates. They must come from stacks the frozen project does not
  hold, so no image enters the union twice, and must name the same image in
  `ptcl2D` and `ptcl3D`.
- A frozen project longer than the current one is accepted when no frozen
  particle lies past the current project's last row.
- A frozen particle rejected in the current project's `ptcl2D` is retired:
  its row stays in place, but it contributes to neither the frozen term nor
  the cohort and remains inactive in the final union. Other frozen particles
  must be active in the frozen project's `ptcl2D`, with the same CTF
  parameters and optics group in both projects.

Every refusal names the first offending particle and is raised before any
project is written. A permutation of rows (a pruned or reordered project) is
refused: `selection` keeps deselected rows by default (`prune=no`), and a
stream must append rows, never renumber them.

Membership is defined once:

```text
frozen = frozen project ptcl3D state > 0 .and. updatecnt > 0
          .and. current project ptcl2D state > 0
cohort = current project ptcl2D state > 0 .and. .not. frozen
```

Frozen-project rows that were selected but never updated (a sampled base run)
join the cohort. The cohort needs at least 5 particles per inherited state
(balanced labelling gives each state `ncohort/nstates`); below 5 % of the
frozen population the run warns and proceeds. An empty cohort and an empty
inherited state are refused.

## 5. Run Structure

- **Masking.** The commander saves the working copy's `ptcl2D` states and sets
  `state=0` in `ptcl2D` and `ptcl3D` for every frozen or retired row, so every
  counting, sampling and labelling routine of the established workflow sees
  the cohort alone, with no add-on branch. The working copy's inherited sigma2
  registration is dropped; the run estimates the cohort's own.
- **Frozen sets.** One frozen accumulation per distinct stage box of the
  inherited ladder, plus the native box, each a `reconstruct3D` on the frozen
  copy with the frozen run's committed sigma2 state (a cropped box uses a
  prefix of its shells). A set is written as `frozen_stateNN_boxBBBB_*` with
  its manifest, under a run context (`solve3D_addon_frozen_context.txt`)
  that records the backend, the state layout, the frozen counts and the row
  counts of both projects. Sets are never clipped or padded.
- **Entry.** Stage 3 (`PROB_REFINE_STAGE`) with `pgrp_start=pgrp`: no symmetry
  search. The cohort gets random orientations and uniform random labels across
  every inherited state; its `res` and `res05` are cleared. The native-box
  frozen-only maps are the stage-3 references, trusted as they are (no CC
  pose initialisation); their correlation with the frozen run's final maps is
  logged as a provenance check and warned about below 0.9.
- **Ladder and limits.** The planned ladder (limit and crop per stage) is the
  frozen run's, from its manifest, up to its last stage; the add-on cannot
  re-plan, lower or raise it. The stage limits follow the `solve3D` rule:
  FSC=0.5 promotion at every stage boundary past `FSC05_PROMOTE_MIN_STAGE`
  from the FSC measured on the union, and the NU handoff in the NU stages.
  Each stage boundary logs the add-on's limits next to the frozen run's.
- **Early stopping.** The add-on context switches stage-3 early stopping on
  (`overlap`, default 0.95); the stage controller is otherwise unchanged.

## 6. Sampling, Trailing and Multi-State

- **Sampling** runs unchanged on the masked cohort: `nptcls_eff`, the
  full-sampling switch (`nsample/cohort > 0.9`), `update_frac` and the
  realized fractions are the cohort's. No frozen index enters a sample, a
  probability table or a realized fraction.
- **Trailing** applies to the cohort alone. Under a sampled cohort, per state
  `s` and half `h`:

  ```text
  T_C(s,h,t) = (u_s/f_s) * P_C(s,h,t) + (1 - u_s) * T_C(s,h,t-1)
  U(s,h,t)   = F(s,h) + T_C(s,h,t)
  ```

  The cohort chain is written before the frozen term `F` is added; restoration,
  FSC, priors and NU filtering consume `U`. Before the first trailing stage the
  stage-boundary reconstruction seeds a full-mass, cohort-only chain through
  `trail_seed`; the frozen term never enters the chain.
- **Multi-state.** One frozen state gives `single`, more gives `independent`,
  whatever the frozen run's `multivol_mode` was (`base_multivol_mode` and
  `split_stage` are provenance only). Frozen rows keep their labels; the
  `independent` `prob`/`prob_neigh` policies update the cohort; no consensus
  accumulator, split, `prob_state` or docked neighbourhood is built. Every
  inherited state is carried: a state without cohort particles is
  reconstructed from its frozen term.

## 7. The Frozen Term

The frozen term is summed into every stage reconstruction, with coefficient
one, before any restoration or prior:

```text
gridding:  S_eo = S_eo(cohort) + S_eo(frozen),  rho_eo = rho_eo(cohort) + rho_eo(frozen)
pcg:       B = B(cohort) + B(frozen),           D = D(cohort) + D(frozen)   (before end_accum)
```

Gridding adds it in `restore_state_from_parts`, after the trailing blend; PCG
in the master's half job, in both execution modes, after the chain write. It is added before every zero-current early-out, so a state or half
with no cohort contribution is assembled from the frozen term alone. Consumers
open frozen sets only through the internal `frozen_rec` handshake, which names
the run context and cannot be set on a command line; every set is validated
against the context and the consumer's grid before a byte is read, and a set
of another reconstruction weighting (euclid or cc) is refused.

The frozen term serves the stage references only. Its sigma2 weighting is the
frozen run's, the cohort's is the add-on's, so the union during the stages is
a cohort-specific weighting model; the final map is not (section 8).

## 8. Final Reconstruction and Chaining

Before its final reconstruction the add-on restores the frozen rows (the
frozen project's 3D records through `transfer_3Dparams` plus the state, and
every saved `ptcl2D` state) and drops the cohort-only sigma2 registration. The
consumability check does not look at the active set, so a cohort-only state
would otherwise weight the union with the cohort's noise model.

`calc_final_rec` then reads every particle, finds no consumable sigma2 state
and bootstraps the union's exactly as `bootstrap_rec3D` does for any project
without consumable sigmas: the image-power seed of every particle, the
gridding ML bootstrap map, one residual pass that commits the canonical state
at native sampling, and the shipped map on it. The residual pass leaves the
particle field as it found it ([refine3D_policy.md](refine3D_policy.md),
section 5.1), so the frozen rows stay exactly the frozen project's.

The output therefore carries the union's committed residual sigma2 state, and
its manifest is eligible: it is the frozen input of the next add-on. Add-ons
chain, each frozen on the previous output and searching only the particles
appended since; a stream never rebases for the sake of its sigmas. The next
add-on's frozen accumulations at cropped stage boxes use a prefix of the
state's shells.

## 9. Distributed Execution

Distributed runs (`nparts > 1`) split the rows into contiguous partitions
that balance the particles with state > 0, as every distributed 3D run does
(`split_nobjs_active`, through `qsys_env%new(..., l_active)` in refine3D,
reconstruct3D, `prob_align` and `prob_align_neigh`). Partition *k* ends on the
last active row of the *k*-th share of the even split of the active rows, so
the masked frozen rows join the partition of the next active row and each
partition holds an equal share of the cohort. A stream's frozen rows are its
first rows: the first partition spans all of them plus its share (HolJunk
update 2: rows 1-53652, then 4127 new rows per partition). An interleaved
selection gets boundaries that follow its active rows. With every row active,
as on the restored union of the final reconstruction, the partitions are the
even split.

No consumer assumes a particular split: the master merges the orientation
documents by each document's own range (`merge_algndocs`), and the sigma2,
probability-table and power-spectrum merges read the ranges from their files.

A partition without sampled particles remains a valid transaction: with fewer
active rows than partitions, or when a sample drawn over the whole project
(class-balanced or probabilistic) misses a partition. In a distributed worker
(`part` set) every per-partition sampler returns an empty sample instead of
stopping (`allow_empty`: probabilistic reproduction, class-balanced,
update-count and fill-in sampling); a shared-memory run still stops on an
empty sample. `prob_tab` and `prob_tab_neigh` write a table without
candidates, which `prob_align` reads; the refine3D matcher emits the unchanged
committed sigma2 slice, the range's orientations, zero PCG accumulators when
partial reconstructions are written, and `JOB_FINISHED`.

The local queue checks no exit status: a worker that stops without
`JOB_FINISHED` leaves the master waiting.

## 10. Outputs and Publication

- `mkdir=yes` (default): the run works on a copy of the current project in a
  new run directory. When the run has completed, the finished project replaces
  the current project file (written beside it, then renamed over it) and
  registers the run manifest by absolute path; a failed run never touches the
  current project.
- `mkdir=no` (NICE, the stream): the current project is the working project
  and the run directory is the project's directory.
- The run directory holds the union maps and final products as `solve3D`
  writes them, the union sigma2 state, `solve3D_manifest.txt` (the add-on's
  own, eligible), `solve3D_addon_report.txt`, the frozen copy in `frozen/`,
  the frozen sets and their run context, and with `addon_diag=yes` the
  cohort-only reconstruction in `addon_diag/`.

## 11. Validation Report

Every run ends with a report against the frozen solution, logged and written
to `solve3D_addon_report.txt`:

- per state, the union and frozen populations;
- per state, the union FSC against the frozen run's: FSC=0.5 and FSC=0.143,
  the move of the FSC=0.143 shell, the mean FSC gain up to the frozen run's
  FSC=0.143 shell, and a verdict: IMPROVED or REGRESSED when the shell moved
  by more than one shell, UNCHANGED otherwise;
- per state, the union map's correlation with the frozen map inside the mask
  up to the frozen run's FSC=0.143 resolution; below 0.9 the union map is
  also docked onto the frozen map, which tells a moved frame from a changed
  structure;
- with `addon_diag=yes`, per state, the cohort-only map (particles the frozen
  run never saw, reconstructed without the frozen term on a copy with the
  frozen rows masked) against the frozen map: its FSC and correlation
  cross-validate the cohort's alignment;
- both runs' limits for every stage.

A regression is warned about and the result is published all the same;
rejecting it is the user's call.

## 12. Streaming Integration

For `commander_stream_p07_solve3D_multistate` or any driver that grows a pool:

- **Pool per update.** Build each update's current project as the previous
  solution's project (all its rows, in order) with the new classified sets
  appended as rows: their `ptcl2D` records and 2D class labels kept. p07 does
  this by merging each pool publication into its rows and taking the
  publication's class table (`doc/policies/stream/stream_3D_ingestion_policy.md`).
  `merge_projects` is not a substitute: it drops the 2D classification, and
  `solve3D` requires it.
- **Separate files.** The frozen and current projects are different files;
  the add-on refuses aliases. A pool that `solve3D` ran on in place is the
  frozen project; the next update's pool is a new project file.
- **Chaining.** Update 1 is frozen on the base run, update `n+1` on update
  `n`'s output: each update searches only its new sets.
- **Selection.** A later 2D rejection may retire a frozen particle by setting
  its current `ptcl2D` state to zero. Rows must remain present and ordered;
  pruning or renumbering is still refused by identity validation.
- **Execution.** `mkdir=no` in each update's own directory; `nparts`,
  `worker_server` and `worker_priority` pass through.

## 13. Tests

- Unit sub-suites: `project superset` (a 20-row current project with a 14-row
  frozen project, every identity refusal, membership, masking and restoration),
  `solve3D manifest` and `solve3D addon report` in `unit_project`;
  `frozen accumulator` and `volume pair metrics` in `unit_reconstruction`;
  `addon report docking` in `lib_reconstruction`.
- Workflow gate `solve3D_addon` (CTest, highlevel): simulated particles of
  a symmetry-broken 6VXX map in two stacks. The base runs on a 2000-row
  frozen project (the first set, a seeded 75 % selection), the add-on on a
  3000-row current project that appends the second set. It checks the frozen
  inputs unchanged, the frozen rows restored exactly, the union sigma2 state
  registered, the output valid as a frozen input, the report, cohort coverage,
  cohort poses against the truth and the union map against the truth and the
  base map.

## 14. Code Map

| File | Owns |
| --- | --- |
| `src/main/ui/simple/simple_ui_solve3D.f90` (`new_solve3D_addon`), `src/main/ui/simple_ui_program.f90` (`accepts`, `UI_ENVIRONMENT_KEYS`) | the command-line contract |
| `src/main/exec/simple_exec_solve3D.f90` | the program entry |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | `commander_solve3D_addon` (key check against the UI entry, manifest and identity validation before any write, the command line for `exec_solve3D`, publication) and the add-on route of `exec_solve3D` (prologue, frozen sets, stage entry, restore, epilogue, manifest) |
| `src/main/solve/simple_solve3D_manifest.f90` | the run manifest: writing, reading, `validate_frozen`, `replay`, the recorded inputs (`manifest_records_input`) |
| `src/main/project/simple_project_superset.f90` | identity, membership, masking and restoration |
| `src/main/volume/simple_frozen_accum.f90` | the frozen sets and their run context |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90`, `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | the frozen add on the gridding and PCG backends |
| `src/main/solve/simple_solve3D_controller.f90`, `src/main/solve/simple_solve3D_utils.f90` | the add-on context (stage-3 early stopping, the `frozen_rec` handshake), `calc_frozen_rec`, the ladder from the manifest |
| `src/main/solve/simple_solve3D_addon_report.f90`, `src/main/volume/simple_volpair_metrics.f90` | the validation report |
| `src/main/simple_final_rec.f90`, `src/main/commanders/simple/simple_commanders_refine3D.f90` (`exec_bootstrap_rec3D`) | the final reconstruction that bootstraps the union's sigma2 state |
| `src/utils/simple_map_reduce.f90` (`split_nobjs_active`), `src/utils/qsys/simple_qsys_env.f90` (`new`, `l_active`), `src/main/project/simple_sp_project_core.f90` (`merge_algndocs`) | partitions that balance the active particles, and the range-checked merge |
