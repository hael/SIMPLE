# Sequential per-state refinement of multi-state projects (state sweep)

Status, 2026-11-27: planned. Base: `fe6a07c22`. Supersedes: none. Depends on: the `state=X` work-project
path of `refine3D_auto` and `refine3D_local` (validated: yes, `flex_pca_blobs` runs it per state, hard and
weighted); the frozen state weight set of `flex_pca` (validated: yes, same entry). Line numbers are at the base.

## 1. Why

A multi-state project is refined today one state at a time, by hand: `refine3D_auto|local state=X` builds a
single-state work project, refines it and leaves the result there; the parent is never written
(`refine3D_policy.md`, "Per-state M-estimation"). Nothing folds the states back into one project, so a
follow-up run (`refine3D_local` after `refine3D_auto`, a second local round) cannot start from the parent,
and under `m_estimator=flex` a particle refined in several states has one pose per state with nowhere to keep
them. The joint alternative (one multi-state refine3D whose states do not interact) was weighed on 2026-11-27
and dropped: running each state "like a single-state refinement" needs per-state sigma2 spectra and matching
bands throughout the polar, Cartesian, sigma2-store and reconstruction code (about 15-20 files). The
single-state path already has both. This plan runs it once per state and merges the results.

## 2. Contract

- **Off-state invariant.** With `state=` given, or a project with at most one populated state (hard) /
  one state of applied mass (flex), `exec_refine3D_workflow` takes the existing path; no sweep code runs
  beyond the one predicate that decides. The `state=X` path stays byte-identical at one thread (Phase 2 gate).
- **Children are ordinary runs.** Each state is refined by the unchanged `state=X` path, in-process, in its own
  numbered directory below the sweep's directory. The sweep adds no scientific behaviour: no shared sigma2,
  band, sample, convergence or reference between states.
- **One project owner.** The sweep owns its directory's copy of the parent and is the only writer of it; each
  child owns its work project. The project given on the command line is never written (as `state=` today and
  `solve3D_chunks`).
- **Labels are the parent's.** The sweep never changes a `ptcl3D` state label. Under `m_estimator=flex` the
  weights are frozen, as today.
- **Merged parent is a valid source.** After a completed sweep, `refine3D_local` (sweep or `state=X`) starts from
  the parent without a provenance refusal, and under flex every state starts from its own refined poses.

## 3. Decisions and rulings

1. (2026-11-27, Hans) **Sequential sweep, not a joint engine.** States run one after another; each child gets
   the full `nthr`/`nparts` of the command line. No concurrency between states.
2. (2026-11-27, Hans) **`m_estimator=flex` is supported by the sweep**, in `refine3D_auto` and `refine3D_local`.
   `refine3D_states` keeps refusing it (a search over states moves the labels the weights were computed for);
   its refusal message points to the sweep instead of `refine3D_auto state=X`.
3. **Trigger.** `refine3D_auto|local` without `state=` on a project with more than one state to refine sweeps.
   States to refine: hard, every `ptcl3D` label with population > 0; flex, every state of the weight set with
   applied mass > 0 (`get_applied_weights(s, 0.)`). Ascending order. No new key; no subset selection.
   `vol1` with a sweep is refused (it is one state's map). With one such state the existing path runs as now
   (decision 3 is the off-state predicate).
4. **Placement.** The branch sits in `exec_refine3D_workflow` after the defaults block (lines 200-243) and before
   `controller%new` (line 244), so the sweep's `params` carries the same settings as its children
   (`settings_identity` of the provenance records, decision 10).
5. **Re-entry.** The sweep is a controller type `refine3D_state_sweep` (`src/main/simple_refine3D_state_sweep.f90`,
   `new`/`run`/`kill`). The commander hands itself to it: `sweep%run(params, cline, self)` with
   `class(commander_base)`, and the sweep calls `cmd%execute(cline_state)` with `state=s mkdir=yes`. No module
   cycle. `exec_refine3D_auto`, `exec_refine3D_local` and `exec_refine3D_workflow` become `recursive` (the build
   has `-fopenmp`, which implies `-frecursive`; the attribute states it). The child cline is a copy of the
   entry cline after the defaults block plus `state` and `mkdir=yes`. After each child the sweep returns to its
   directory with `simple_chdir` and resets `CWD_GLOB` (`simple_chdir` does not, `simple_syslib.f90:869`).
6. **Hard merge.** Work row `i` goes to parent row `parent_rows(i)`: the whole `ptcl3D` record is copied, then
   the parent's `state` label is restored. The parent-to-work row map and its stack/index check, today inline in
   `prepare_refine3D_state_work_project` (lines 470-481), move to one `sp_project` routine used by the
   preparation and the merge. The row selection stays `select_refine3D_state_rows` (`simple_refine3D_strategy.f90:202`).
7. **Flex merge.** For state `s`: rows whose parent label is `s` take the whole record as in decision 6; every
   selected row's pose (e1, e2, e3, x, y) goes to column `s` of the state pose set (decision 8). A particle's
   `ptcl3D` pose is therefore always its hard-labelled state's pose.
8. **State pose set.** New class `state_pose_set` (`src/main/project/simple_state_pose_set.f90`), file layer
   `src/fileio/simple_state_poses_file.f90`, built like `state_weight_set`: one file per state of `nptcls` rows
   of five `real32` (e1, e2, e3, x, y) and one `int8` flag (refined for this state), a manifest published last,
   an out-segment entry `imgkind=state_poses` pointing at the manifest. Identity: particle layout digest,
   the weight set's generation and layout digest, `nstates`. Written only under flex.
9. **Pose-set validity and seeding.** On open the set must match the current weight set's identity and agree
   with `ptcl3D` at every flagged row whose label is the column's state (exact `real32` equality: the values are
   copies). Any mismatch: one warning, the set is ignored, seeds come from `ptcl3D` (the state after `flex_pca`).
   One routine seeds a work project's rows from column `X`; both `prepare_refine3D_state_work_project` and the
   local preflight (`preflight_cont_startup`, the `check_cont_seeds` of the selected rows) call it, so a
   `state=X` run outside a sweep seeds the same way. A hard sweep removes a registered pose set from the parent.
10. **Out segment and provenance.** The work project's state-1 entries (`vol`, `fsc`, `vol_cavg` when present)
    are registered as state `s` with their fields (`pop`, `mass`, `ess`); paths are already absolute. A new
    `sp_project` routine beside `remove_state_artifacts_from_osout` does it. After the last state, each state's
    reference record is published again against the merged parent (`publish_reference_set`, in place beside the
    child's final map; map, halves, FSC and `align_lp` from the child's own record). A child without a valid
    record gets none; a later local run then takes its documented legacy-pair route. The superseded work projects
    no longer validate against their records; continuation is from the sweep's project.
11. **Write cadence.** The parent's `ptcl3D` and `out` are written after each state. Pose-set state files are
    written per state, the manifest after the last state, the provenance records after the manifest. An
    interrupted sweep leaves a consistent parent with the finished states merged and no pose set.
12. **Not carried back:** trailing chains (a later local run seeds a new chain with a full first update, the
    documented fallback), sigma2 stores (the parent's store is left; consumers validate it and re-bootstrap on a
    mismatch), sampling bookkeeping across states. The single `vol_msk` entry of the out segment
    (`add_vol2os_out`) keeps the last state's mask; each state's masks stay in its reference record.
13. (2026-11-27, Hans) **`refine3D_local refine=prob_neigh` requires `cont_view`.** An unset `cont_view` becomes
    `lm`; `cont_view=no` is refused. The local `nspace=20000`/`nspace_sub=500` defaults (lines 149-150) go; the
    `cont_view` default of `nspace=1000` stays. This is the maintainer's switch of `polar_continuous_view.md`
    decision 6 (pending the PorF read in its report, section 8); the engine's `prob_neigh` without `cont_view`
    stays (base `refine3D`, `classify3D_refs`).

Defaults taken without a question: the type and file names above; the sweep's log line per state (state, rows,
mass and ESS under flex, run directory, seconds) and a `refine3D_state_sweep.txt` table (state, run directory,
work project) in the sweep's directory; ascending state order.

## 4. Not in this plan

- A joint multi-state engine, per-state sigma2 spectra or per-state matching bands.
- Concurrency between states; restart or `continue` of an interrupted sweep (rerun from the parent).
- Carrying trailing chains or sigma2 stores back (decision 12).
- Any change to `refine3D_states`, `refine3D_states_cavgs`, `flex_pca` or `solve3D` beyond decision 2's message.
- `src/main/strategies/search/**`, `src/main/pftc/**`, `src/main/cftc/**`, `src/main/sigma2/**`: out of bounds.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 1 | result-changing (local `prob_neigh` default) | fast gate; `unit_core`, `unit_cart_align3D`; the new `cont_view` case (Phase 1 exit) |
| 2 | identical | fast gate; `unit_project`; byte identity of the direct `state=X` runs (section 6) |
| 3 | new path, off-state identical | Phase 2's identity; sweep child ≡ direct run; `unit_project` |
| 4 | new path (flex), off-state identical | as 3, plus `unit_project` pose-set tests |
| 5 | new path (provenance, chain) | the chained local gate (section 6); `unit_reconstruction` |
| 6 | docs, report | fast gate |

## 6. Acceptance

- **Fixture.** `flex_pca_blobs` (`simple_commanders_test_highlevel_flex.f90`): two blob states, a frozen weight
  set, and `run_state_refinements` already refining each state directly with `m_estimator=flex` and `no`
  (`STATE_REFINE_MAXITS = 5`, gate `STATE_REFINE_MARGIN = 0.01`, measured on the profile of 2026-10-19).
- **Determinism** (measured, Phase 0): two direct `refine3D_auto state=X` runs per mode at `nthr=1`, fixed seed.
  If final maps and work-project `ptcl3D` are byte-identical, identity gates below are byte gates. If not, the
  tolerance is twice the largest difference between the two base runs, frozen in Phase 0's Progress entry.
- **Phase 2 identity.** The direct `state=X` runs on the Phase 2 build equal the base runs (maps, `ptcl3D`
  records). Exclusions: paths, `projinfo`, timestamps, logs.
- **Sweep child ≡ direct run** (Phases 3-4). Each child's final map and work `ptcl3D` equal the direct run of
  the same state and mode, under the same rule and exclusions.
- **Merge.** In the swept parent, every work row equals its parent row (whole record, label excepted) for
  hard; for flex, label-`s` rows equal state `s`'s work rows and the pose-set column `s` equals every selected
  work row's pose. Out entries of state `s` name the child's map and FSC with the child's `pop`/`mass`/`ess`.
- **Chain** (Phase 5). Base reference, measured in Phase 0: direct `refine3D_auto state=X`, then
  `refine3D_local` (`refine=cont`, `maxits=3`) on its work project, per state and mode, twice; floor = base
  median map-to-truth correlation minus twice the two-run spread (measured, frozen in Phase 0). Gate: an auto
  sweep followed by a local sweep on the merged parent starts without a provenance refusal and every state's
  local map meets the floor.
- **Reported, not gated:** the per-state FSC curves of the chained run (`fsc_state01.bin` of each child) and the
  sweep table, collected under the entry's run directory.

## 7. Common rules and stop conditions

- No new command-line parameter. No floor, tolerance or assertion relaxed. Fix as you go, with a file-table row
  first. Fast gate runs at one thread; the high-level entry at its registered thread count except the identity
  runs (`nthr=1`).
- **Stop** if: a child changes a `ptcl3D` state label (decision 3's premise; Phase 0 checks it); Phase 0's
  direct runs are not reproducible at one thread and the two-run spread exceeds `STATE_REFINE_MARGIN`; the
  provenance republish needs a change to the record format; an edit is needed outside the file table or in an
  out-of-bounds folder; `flex_pca_blobs` would exceed its registered budget (7,200 s) by more than the measured
  sweep time.
- **No ruling waives:** the fast gate; Phase 2's byte identity; the sweep-child identity; the chained local gate.

## 8. Phases

- **Phase 0, baseline.** Build Release (and Debug for the bounds check), fast gate, record machine and load.
  Run `flex_pca_blobs` twice. Run the determinism and chain references of section 6 on the base. Check that no
  child changes a label (diff parent and work-project labels of the selected rows). List the `ptcl3D` keys a
  child writes. Exit: Progress holds the frozen tolerance or the byte verdict, the chain floor per state and
  mode, the label check, the key list.
- **Phase 1, local `prob_neigh` requires `cont_view`** (decision 13). Commander branch, UI help, parameters
  tester, `refine3D_local` tester. Add a `refine=prob_neigh` (`cont_view=lm` by default) local case to
  `cont_refine3D_1jxy` (none runs `cont_view`, `polar_continuous_view_report.md` section 7).
  Exit: `cont_view=no` with `refine=prob_neigh` refused (tester case); the new case passes; `rg "20000"`
  finds no local default in `simple_commanders_refine3D.f90` outside `refine3D_auto`.
- **Phase 2, row map extracted** (decision 6). The parent-to-work map and its check become one `sp_project`
  routine; `prepare_refine3D_state_work_project` calls it. Exit: `unit_project` test of the map (a reordered or
  foreign stack is refused); Phase 2 identity.
- **Phase 3, hard sweep** (decisions 3-6, 10 first half, 11, 12). The sweep type, the commander branch, the
  absorb and out-segment routines. Extend `run_state_refinements` with an `m_estimator=no` sweep of the parent.
  Exit: sweep-child identity and the hard merge checks; the off-state predicate unit-tested (one populated state,
  `state=` given); `refine3D_state_sweep.txt` lists both states.
- **Phase 4, flex sweep and pose set** (decisions 2, 7-9). The pose set and its file layer; seeding in work-project
  preparation and local preflight; the hard sweep's removal of a registered set. Exit: `unit_project` pose-set
  tests (round trip, identity mismatch ignored, a disagreeing flagged row ignored, manifest-last publication);
  a flex sweep in `run_state_refinements` meets the sweep-child identity and the flex merge checks.
- **Phase 5, provenance and the chain** (decision 10 second half). Republish per state after the merge; add the
  chained auto-sweep → local-sweep case. Exit: the chain gate; `unit_reconstruction` passes.
- **Phase 6, docs and report.** Living descriptions: `refine3D_policy.md` ("Per-state M-estimation", new "State
  sweep"), `refine3D_local_policy.md`, `refine3D_auto_policy.md`, `heterogeneity/refine3D_states_policy.md`
  (decision 2's message), `how2_process_heterogeneous_datasets.md`. Report `refine3D_state_sweep_report.md`.
  Exit: fast gate; every policy above names the sweep; the report exists.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | branch, `recursive`, `prob_neigh` defaults, row-map call, seeding call | 1-4 |
| `src/main/simple_refine3D_state_sweep.f90` (new) | the sweep controller | 3-5 |
| `src/main/simple_refine3D_workflow_controller.f90` | local preflight seeds from the pose set | 4 |
| `src/main/project/simple_sp_project.f90`, `simple_sp_project_ptcl.f90`, `simple_sp_project_out.f90` | row map, absorb, out-entry relabel | 2-4 |
| `src/main/project/simple_state_pose_set.f90` (new), `simple_state_pose_set_tester.f90` (new) | pose set | 4 |
| `src/fileio/simple_state_poses_file.f90` (new) | pose-set files and manifest | 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | `select_refine3D_state_rows`, if needed | 2-4 |
| `src/main/volume/simple_refine3D_provenance.f90` | reading a child's record for republication, if needed | 5 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` (`exec_refine3D_states` message) | decision 2 | 4 |
| `src/main/params/simple_parameters_phases.f90`, `simple_parameters_tester.f90` | `cont_view` cases | 1 |
| `src/main/ui/simple/simple_ui_refine3D.f90` | help of `refine`, `cont_view`, `state`, `m_estimator` | 1, 3 |
| `src/main/commanders/simple/simple_refine3D_local_tester.f90` | local `prob_neigh` cases | 1 |
| `src/main/commanders/test/simple_commanders_test_highlevel_flex.f90` | sweep, merge and chain cases | 3-5 |
| `src/main/commanders/test/simple_commanders_test_highlevel_cont.f90` | local `cont_view` case | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | pose-set suite row | 4 |
| `src/main/commanders/test/simple_test_registry.f90`, `production/CMakeLists.txt` | `flex_pca_blobs` budget, if needed | 3-5 |
| `src/CMakeLists.txt` | new sources | 3, 4 |
| `doc/policies/3D/refine3D_policy.md`, `refine3D_local_policy.md`, `refine3D_auto_policy.md`, `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/how2s/how2_process_heterogeneous_datasets.md` | living descriptions | 1, 6 |
| `doc/implementation_notes/planned/refine3D_state_sweep*.md` | this plan, Progress, report | all |

## Progress

(one entry per phase, written by the run)
