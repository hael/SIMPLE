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
maps in one stage, the last stage the frozen run ran (recorded in its
manifest), preceded by a registration pass. The result is one ordinary project
with every particle posed.

Two uses are supported:

- **A more permissive selection of a fixed project.** The frozen solution was
  computed on a harsh selection; the add-on adds what a more permissive
  selection of the same particles contains. The two projects have the same
  rows and differ in their selection.
- **A growing set in an offline chain.** Each add-on of `solve3D_chunks`
  appends new particles to the previous step's and is frozen on the previous
  step's output, so it searches only the new particles (section 12). The stream
  no longer runs add-ons (13 October 2026): its 3D stage assigns class averages
  to fixed reference volumes.

It is not a refinement of the frozen particles, not a class-average route
(the add-on aligns particles only and never runs `solve3D_cavgs`) and not
an ini3D phase.

## 2. Command Line

The add-on runs with the frozen run's settings; its command line carries
what governs compute, convergence and diagnostics, plus three settings of the
frozen run it may override because they describe the particles being added
rather than the solution. The contract is the
program's UI entry (`new_solve3D_addon` in `simple_ui_solve3D.f90`),
checked by the wrapper through `ui_program%accepts`
([ui_layer_policy.md](../ui_layer_policy.md)); there is no second key list:

| Group | Inputs the UI entry declares |
| --- | --- |
| Projects | `projfile_frozen` (`projfile` and `mkdir` as for every project program) |
| Compute | `nparts`, `nthr` |
| Sampling and convergence | `nsample` (default: the frozen run's effective value), `overlap` (default 0.95, applies only when the last stage is stage 3) |
| Overridable settings of the frozen run | `balance`, `nclust`, `mskdiam` (default: the frozen run's value, replayed from its manifest) |
| Diagnostics | `addon_diag` |

The execution environment every program accepts from its launcher passes
through unchanged (`UI_ENVIRONMENT_KEYS`: the queue system, NICE, and the
stream's `worker_server` and `worker_priority`). Every accepted key is
forwarded as given.

Refusal is by key, so an inherited value cannot be confirmed silently. A key
the frozen run's manifest records as one of its inputs
(`manifest_records_input`: the solution, reconstruction and search policy, and
the entry routes as provenance) is refused as set by the frozen run, except
the overridable keys (`manifest_overridable_input`: `balance`, `nclust`,
`mskdiam`), which the UI entry accepts: the replayed value is their default,
a value on the add-on command line replaces it, and the wrapper logs the
replacement against the frozen run's value. They govern how the cohort is
sampled and masked, which depends on the particles being added: the frozen
accumulators are mask-free raw sums, and the mask enters at restoration, in
the solve support and in matching, on the union. The add-on's own manifest
records the values it ran with, so a chained add-on inherits the override.
`nu_stage` is a recorded and replayed input: the add-on starts NU filtering at
the stage its base run did, so its last stage is filtered, or not, as in the base.
Any other key the UI entry does not accept is refused as not an input.

A superset project without a 2D solution (its `ptcl2D` never searched) runs
under `balance=none` only, given on the add-on command line when the frozen
run used `class` or `cavg`; the selection is still its `ptcl2D` state flags
and the ladder is the frozen run's, so nothing else changes
([solve3D_policy.md](solve3D_policy.md), "Without a 2D solution").
Ordinary `solve3D` and `solve3D_cavgs` refuse `projfile_frozen` and
`addon_diag`.

`center=no` is forced: a re-centred reference would map a shift onto cohort
particles that the frozen accumulators never receive.

`rec_backend=gridding` is forced. Add-on chains look for initial solutions and
speed, not the precision of the kernel preconditioned conjugate-gradient (PCG)
reconstruction, so they reconstruct by gridding only (ruled 2026-10-13): the
UI entry offers no PCG controls, the backend is not replayed from the
manifest, and `rec_backend` or a PCG control (`maxits_pcg`, `maxits_ml`,
`pcg_solvent`, `pcg_solvent_lambda`, `pcg_solvent_check`) on the command line
is refused by name.

## 3. The Frozen Input

The run manifest is the only route into the add-on. `exec_solve3D` writes
`solve3D_manifest.txt` at the end of every completed run (after the final
reconstruction) and registers it in `projinfo`;
a single-state run stopped before stage 3 (stage 4 if it searches a symmetry
axis) writes none and is not a frozen input; a multi-state run always writes
one, and the add-on refuses it when its last stage is before stage 3. A run
stopped explicitly at that stage or later
(`nstages`, e.g. a single-state `nstages=5` base) does the final reconstruction
and writes one; the add-on then runs that stage. The frozen project
must be an eligible `solve3D` or `solve3D_addon` output whose manifest
validates against the project that registered it:

- the run identifier registered in `projinfo` matches the manifest's;
- row count, particle layout digest, stack table and optics/CTF parameters
  are the ones the manifest recorded;
- every state map registered in the project is the map the run produced
  (digest);
- the committed residual sigma2 state has the recorded committed header and
  grouped-section identity; manifest validation never scans its particle spectra;
- the current project's native box and sampling equal the frozen solution's;
- the frozen run reached stage 3 or later, and its last stage's box does not
  exceed the native box;
- the frozen run reconstructed by gridding: a base run whose manifest records
  `rec_backend=pcg` is refused, naming its backend, before anything is
  written (`solve3D_manifest%get_rec_backend`; no recorded backend is the
  gridding default). A `solve3D` base run keeps both backends.

The manifest reader keeps no backwards compatibility (release 4): a field
or input key it does not know, such as those of removed features
(`ptcl_src`, `objfun_den`, `conical_fsc`, `inpl_cont`), makes the manifest
unreadable, so a solution from an older build cannot be extended.

The frozen project is never written. The add-on works on a copy of it in the
run directory (`frozen/`, with its sigma2 state), whose sigma2 state must retain
the frozen run's committed header and grouped-section identity through every
frozen accumulation. This is intentionally not a byte-equality check over the
per-particle spectra.

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
  its row stays in place and remains inactive in the final union. The first
  add-on of a chain leaves it out of the frozen term it accumulates; once a
  particle is inside carried sums (section 8) it stays in them, an accepted
  consequence (ruled 2026-10-13): it is inactive in the project but stays in
  the frozen volumetric components already generated. A row inside the frozen
  run's carried sums that is not frozen (retired earlier) and that the current
  project selects again is retired too: searching it would count it twice.
  Other frozen particles
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
  the cohort alone, with no add-on branch.
- **Sigma2, copy-on-write.** The inherited registration is dropped, and the
  run opens a sigma2 state of its own in its directory
  (`sigma2_state_extend`): the frozen run's committed records of the frozen
  and retired rows are copied unchanged, every cohort row starts from the
  frozen global curve of its half (the project's half-set partition is made
  first), and both global means are taken over the cohort, so they start as
  the frozen means. Nothing is seeded from image power. The stage updates the
  cohort's records only; the means it commits are cohort means while the
  frozen rows are masked. The inherited state belongs to the frozen output and
  is never written (rollback).
- **Frozen sets.** One set per state at the last stage's box, plus the native
  box, under a run context (`solve3D_addon_frozen_context.txt`) that records
  the backend (always gridding), the state layout, the frozen counts and the
  row counts of both projects. A set is written as `frozen_stateNN_boxBBBB_*`
  with its manifest. Sets are never clipped or padded. Their source depends on
  the frozen run (section 8):
  - *The first add-on of a chain* (the frozen run is a `solve3D`, or records
    no carried bundle) reads the frozen particles once per box: a
    `reconstruct3D` on the frozen copy with the frozen run's committed sigma2
    state (a cropped box uses a prefix of its shells).
  - *Later add-ons* adopt the frozen run's carried sums: each carried set is
    validated against the bundle (headers and sizes) and copied under the run
    context; no particle is read. The counts are the bundle's, retired
    particles included.
- **Entry.** `pgrp_start=pgrp`: no symmetry search. The cohort gets random
  orientations and uniform random labels across every inherited state; its
  `res` and `res05` are cleared. The native-box frozen-only maps are the
  references, trusted as they are (no CC pose initialisation). On the first
  route they are restored from the frozen sets and their correlation with the
  frozen run's final maps is logged as a provenance check (warned below 0.9);
  on the carried route they are the frozen run's final maps and halves, which
  are the restoration of the carried native sums.
- **Stage policy.** The add-on runs two refine3D steps, nothing before them
  and nothing between:
  1. *Registration pass.* One full-cohort iteration of global `refine=prob` at
     `nspace=1000` (`NSPACE_REGISTER`), at the last stage's box and planned
     limit, with the shift range of the first prob stage. It is unfiltered and
     unmasked, does not trail, and updates every cohort particle, so the frozen
     term enters at full mass. It brings the random cohort poses into register
     with the frozen references; a cohort particle with no previous alignment
     needs no more than this to be found by the stage's own search.
  2. *The last stage the frozen run ran*, with the controller's configuration
     for that stage unchanged: its search mode (`prob_neigh` is a global pass
     over `nspace_sub` directions, so the registered poses are a start, not a
     constraint), `nspace`, iteration budget, early stopping, sampling, trailing,
     NU filtering and masking. The trailing chain, where the stage trails, is
     seeded from the registered poses at the stage boundary.
- **Ladder and limits.** The planned ladder (limit and crop per stage) is the
  frozen run's, from its manifest; the add-on cannot re-plan, lower or raise
  it, and uses its last stage only. The stage limit follows the `solve3D`
  rule: FSC=0.5 promotion at the stage boundary from the FSC measured on the
  union, and the NU handoff when the stage is an NU stage. The registration
  pass keeps the planned limit (no promotion). The stage logs the add-on's
  limits next to the frozen run's.
- **Early stopping.** The add-on context switches stage-3 early stopping on
  (`overlap`, default 0.95) when the last stage is stage 3; the stage
  controller is otherwise unchanged.
- **Manifest.** The stage is read from the frozen run's manifest
  (`last_stage`, with the whole planned ladder). The add-on's own manifest
  records that stage as its first and last stage, so a further add-on on the
  result runs the same stage.

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
  FSC, priors and NU filtering consume `U`. Before the stage, if it trails, the
  stage-boundary reconstruction (after the registration pass) seeds a full-mass, cohort-only chain through
  `trail_seed`; the frozen term never enters the chain.
- **Multi-state.** One frozen state gives a single-state add-on, more gives
  the multi-state (independent) policies, decided by the frozen run's
  `nstates` alone. Frozen rows keep their labels; the multi-state
  `prob`/`prob_neigh` policies update the cohort; no consensus accumulator or
  split is built. Every
  inherited state is carried: a state without cohort particles is
  reconstructed from its frozen term.

## 7. The Frozen Term

The frozen term is summed into every stage reconstruction, with coefficient
one, before any restoration or prior:

```text
S_eo = S_eo(cohort) + S_eo(frozen),  rho_eo = rho_eo(cohort) + rho_eo(frozen)
```

Gridding adds it in `restore_state_from_parts`, after the trailing blend; the
frozen sets have no PCG side, and the PCG backend refuses the frozen-set
handshakes. It is added before every zero-current early-out, so a state or half
with no cohort contribution is assembled from the frozen term alone. Consumers
open frozen sets only through the internal `frozen_rec` handshake, which names
the run context and cannot be set on a command line; every set is validated
against the context and the consumer's grid before a byte is read, and a set
of another reconstruction weighting (euclid or cc) is refused.

The frozen term keeps the weighting it was accumulated with: the frozen run's
sigma2 state for the first add-on's own accumulation, each earlier run's for
the carried sums. Contributions are never re-weighted when sigmas change; that
old weighting is an accepted error (ruled 2026-10-12).

## 8. Final Reconstruction and Chaining

The final reconstruction is incremental (`calc_final_rec_addon`): it reads the
run's own particles (the cohort, with the frozen and retired rows masked) and
the frozen sets, never the frozen particles, and has no image-power seed.
Between its passes the union's membership is restored and the sigma2 state is
republished with both global means over the union (`sigma2_state_union_means`:
no record changes, the generation advances), then the frozen rows are masked
again:

1. **Bootstrap map:** the frozen set at the native box plus the cohort at its
   final poses, on the stage's sigmas (euclid, ML-regularized gridding, the
   stage's filtering).
2. **Residual pass:** the cohort against the bootstrap map at native sampling
   (`refine=sigma`); it commits the cohort's new sigmas, leaving the frozen
   records as they are. The particle field is left as it was found
   ([refine3D_policy.md](refine3D_policy.md), section 5.1).
3. **Carried sums and the shipped map:** the cohort accumulated again on its
   new sigmas, with the frozen set added, at the stage box first when it
   differs from the native box, then at the native box, whose map is shipped.
   Each reconstruction writes the union's raw sums (frozen set plus cohort)
   through the internal `frozen_carry` handshake, before any restoration, as
   the run's carried sets (`carried_stateNN_boxBBBB_*`).

The union is then restored and its maps registered. The pass count is three
when the stage box is the native box and four when it differs; `addon_diag=yes`
adds the cohort-only map.

**The carried bundle.** The carried sets, the rows they contain
(`solve3D_addon_carried_membership.bin`: the state of each row inside the sums,
zero outside) and the bundle manifest (`solve3D_addon_carried_bundle.txt`) form
the run's bundle. The manifest records the run, its parent (the frozen run),
the route (`first` or `carried`), the state layout and row count, the boxes and
their sampling, the weighting, the sigma2 generation of the run's final state
and the particle count per state. It is written last, after every set and the
membership; the solve manifest records it by the digest of this small manifest
(the sums are never digested), and the project registers the solve manifest
only after that. Validation checks headers, sizes, lineage and generation:
another run's, a truncated, a mixed or an incomplete bundle is refused. Each
adopted run owns its complete bundle, so a chain survives the removal of
earlier run directories, and a failed or rolled-back run leaves the previous
bundle untouched.

**Chaining.** The output carries the union's committed sigma2 state (the frozen
records unchanged, the cohort's new ones, means over the union), its eligible
manifest and its bundle: it is the frozen input of the next add-on, which
adopts the bundle and reads only its own particles. Over a chain, each particle
is read a bounded number of times however many add-ons follow. The offline
chains publish their last add-on as it is, with the old weighting of the
frozen contributions; no terminal full reconstruction is added (ruled
2026-10-13).

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
- `mkdir=no` (NICE): the current project is the working project
  and the run directory is the project's directory.
- The run directory holds the union maps and final products as `solve3D`
  writes them, the union sigma2 state, `solve3D_manifest.txt` (the add-on's
  own, eligible), `solve3D_addon_report.txt`, the frozen copy in `frozen/`,
  the frozen sets and their run context, the carried bundle (section 8), and
  with `addon_diag=yes` the cohort-only reconstruction in `addon_diag/`.

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
rejecting it is the caller's call. The stream's multistate stage rolled such a
run back until 13 October 2026, when it stopped running `solve3D_addon`
(`doc/policies/stream/stream_3D_ingestion_policy.md`, section 4, item 11).

## 12. Chains

The stream's multistate stage stopped running add-ons on 13 October 2026; it
now assigns class averages to fixed reference volumes
(`doc/policies/stream/stream_3D_ingestion_policy.md`). The add-on passes serve
the offline chain below. A driver that grows a pool of particles follows these
rules:

- **Pool per update.** Build each update's current project as the previous
  solution's project (all its rows, in order) with the new classified sets
  appended as rows: their `ptcl2D` records and 2D class labels kept.
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

### Automatic balanced chunks

`simple_exec prg=solve3D_chunks projfile=<project> <solve3D inputs>` runs the
same base-then-add-on chain on all selected particles of a project, without
caller-chosen chunk sizes. Every run samples 5000 particles per iteration in
total and searches 10000 to 30000 new particles. Its runs are consecutive
whole stacks in arrival order, every add-on is adopted, and the final union is
published over the workflow's project copy. The policy is
[solve3D_chunks_policy.md](solve3D_chunks_policy.md).

## 13. Tests

- Unit sub-suites: `project superset` (a 20-row current project with a 14-row
  frozen project, every identity refusal, membership, masking and restoration),
  `solve3D manifest` and `solve3D addon report` in `unit_project`;
  `frozen accumulator` and `volume pair metrics` in `unit_reconstruction`;
  `addon report docking` in `lib_reconstruction`; `solve3D chunks plan`
  (the `nsample` rule, the run count, balanced whole-stack partitioning, source
  preparation, step projects and reports) in `unit_project`.
- Unit cases of the carried sums: a two-add-on chain of frozen set, carried
  sums and adopted sums equal to the direct union (sums, densities and the
  restored map) and the bundle refusals (`frozen accumulator`); the
  copy-on-write sigma2 extension and the union means (`sigma2 state` in
  `unit_pftc_align2D3D`); a reselected row inside the sums (`project
  superset`); the bundle artifact (`solve3D manifest`).
- Workflow gate `solve3D_addon` (CTest, highlevel): simulated particles of
  a symmetry-broken 6VXX map in three stacks. The base runs on a 2000-row
  frozen project (the first set, a seeded 75 % selection); a first add-on on a
  2500-row project that appends the second set reads the frozen particles once
  and publishes the first bundle; a second add-on on the 3000-row project that
  appends the third set adopts it (route `carried`, parent the first add-on).
  It checks the frozen inputs unchanged, the frozen rows restored exactly, the
  union sigma2 state registered, the outputs valid as frozen inputs, both
  bundles, the report, cohort coverage, cohort poses against the truth and the
  union map against the truth and the base map.

## 14. Code Map

| File | Owns |
| --- | --- |
| `src/main/ui/simple/simple_ui_solve3D.f90` (`new_solve3D_addon`), `src/main/ui/simple_ui_program.f90` (`accepts`, `UI_ENVIRONMENT_KEYS`) | the command-line contract |
| `src/main/exec/simple_exec_solve3D.f90` | the program entry |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | `commander_solve3D_addon` (key check against the UI entry, manifest and identity validation before any write, the command line for `exec_solve3D`, publication) and the add-on route of `exec_solve3D` (prologue, frozen sets, stage entry, restore, epilogue, manifest) |
| `src/main/solve/simple_solve3D_manifest.f90` | the run manifest: writing, reading, `validate_frozen`, `replay`, the recorded inputs (`manifest_records_input`) |
| `src/main/project/simple_project_superset.f90` | identity, membership, masking and restoration |
| `src/main/volume/simple_frozen_accum.f90` | the frozen sets and their run context; the carried sets, their bundle manifest and membership |
| `src/main/strategies/parallelization/simple_volassemble_strategy.f90` | the frozen add (gridding only) |
| `src/main/solve/simple_solve3D_controller.f90`, `src/main/solve/simple_solve3D_utils.f90` | the add-on context (the registration pass, stage-3 early stopping, the `frozen_rec` handshake), `calc_frozen_rec`, the ladder from the manifest |
| `src/main/solve/simple_solve3D_addon_report.f90`, `src/main/volume/simple_volpair_metrics.f90` | the validation report |
| `src/main/commanders/simple/simple_commanders_solve3D_chunks.f90`, `src/main/solve/simple_solve3D_chunks_plan.f90` | `solve3D_chunks`: the automatic orchestration with its source preparation, cumulative step projects, statistics and report leaf utilities |
| `src/main/simple_final_rec.f90` (`calc_final_rec_addon`) | the incremental final reconstruction and the carried sums |
| `src/main/sigma2/simple_sigma2_state.f90` (`sigma2_state_extend`, `sigma2_state_union_means`) | the copy-on-write sigma2 transaction |
| `src/utils/simple_map_reduce.f90` (`split_nobjs_active`), `src/utils/qsys/simple_qsys_env.f90` (`new`, `l_active`), `src/main/project/simple_sp_project_core.f90` (`merge_algndocs`) | partitions that balance the active particles, and the range-checked merge |
