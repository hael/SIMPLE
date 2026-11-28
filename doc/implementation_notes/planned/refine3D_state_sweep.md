# Sequential per-state refinement of multi-state projects (state sweep)

Status, 2026-11-27: planned; rewritten the same day after Hans's ruling 3, which replaces the revision reviewed in
`refine3D_state_sweep_review.md`. Base: `d110e7c27` (this note is committed on top of it; Phase 0 records the commit
the run gets). Supersedes: none. Depends on: the `state=X` work-project path of `refine3D_auto` and `refine3D_local`
(validated: yes, `flex_pca_blobs` refines each state directly, hard and weighted) and the frozen state weight set of
`flex_pca` (validated: yes). Source line numbers are at the base.

## 1. Why

A multi-state project is refined one state at a time, by hand. `refine3D_auto|local state=X` builds a single-state
work project, refines it and leaves the result there; the parent is never written (`refine3D_policy.md`,
"Per-state M-estimation"). Nothing folds the states back into one project. So a follow-up run, such as
`refine3D_local` after `refine3D_auto`, cannot start from the parent. The sweep runs the existing single-state
workflow once per state and merges the results into one project.

## 2. Contract

- **Off-state invariant.** With `state=` given, or with at most one state to refine (decision 4), `refine3D_auto` and
  `refine3D_local` take the existing path. Not one line of the sweep runs on that path.
- **Children are ordinary runs.** Child `s` is `exec_refine3D_workflow` called with the user's command line plus
  `state=s`. That is the same workflow, defaults, work-project preparation, parallelization (`nthr`, `nparts`) and
  M-estimation as a direct `state=s` run. The sweep adds no scientific behaviour: nothing is shared between states.
- **Merge once, at the end.** No project is merged or written while children run. After the last child returns,
  the sweep reads its copy of the parent, merges every work project into it, and writes it once, atomically.
- **Labels and identity are the parent's.** The sweep never changes a `ptcl3D` state label, a row's stack, its index
  in the stack or its half. Under `m_estimator=flex` the weights stay frozen, as today.
- **The input is never written.** The project given on the command line and its registered files are unchanged.
- **The merged parent is a valid source.** `refine3D_local`, by sweep or with `state=X`, starts from it without a
  provenance refusal.

## 3. Decisions and rulings

1. (2026-11-27, Hans) **Sequential, not a joint engine.** States run one after another. Each child gets the command
   line's full `nthr`/`nparts`, parallelized as a single-state job. There is no concurrency between states.
2. (2026-11-27, Hans) **Both hard assignment and `m_estimator=flex`**, in `refine3D_auto` and `refine3D_local`.
   `refine3D_states` keeps refusing flex; its message points to the sweep instead of `refine3D_auto state=X`.
3. (2026-11-27, Hans) **Simple and built from what exists.** The sweep is a sequential loop over the existing
   `state=X` workflow, followed by one merge after all children have completed. No pose sidecar, no new command-line
   key, no change to the CLI's job-description write and no fault-injection machinery.
4. **Trigger.** `refine3D_auto|local` without `state=` sweeps when the project has more than one state to refine.
   - Hard mode: every positive `ptcl3D` label with a population above 0.
   - Flex mode: every state of the weight set with positive applied weight on some row, as selected by
     `select_refine3D_state_rows` (`simple_refine3D_strategy.f90:202`).

   States run in ascending order. `vol1` with a sweep is refused, because it is one state's map. With one state to
   refine, the existing path runs; that includes the existing refusal of flex without `state=`.
5. **Dispatch in the thin commanders.** `exec_refine3D_auto` and `exec_refine3D_local` (lines 56-68) evaluate
   decision 4 on the command as received. When it holds, they call a private `run_refine3D_state_sweep(cline,
   l_local)` in `simple_commanders_refine3D.f90`; otherwise they call `exec_refine3D_workflow`, unchanged. The sweep
   dispatches before the workflow mutates the command, so no command snapshot is needed. Children call
   `exec_refine3D_workflow` directly, so there is no recursion, no new type and no new module.
6. **The sweep, step by step.**
   1. Resolve the source project: `projfile`, or discovery as in `preflight_cont_startup`. Before saving the
      child template or changing directories, resolve relative input-file arguments (such as `pcg_mskfile`)
      against the entry directory using the existing path helpers. Children retain those absolute input paths.
   2. Keep a child template: a copy of the command after the commander's `autoscale` default (lines 59 and 66)
      and the input-path resolution above. Children bypass the commander, so the template carries that default.
   3. Build `params` once with `mkdir=yes`. This creates the sweep's numbered directory and its copy of the parent,
      and leaves `cline`'s `projfile` pointing at that copy, so the CLI later appends its job entry there, as for
      any program.
   4. For each state `s`, make a fresh copy of the template and set `projfile` to the copy's absolute path,
      `state=s` and `mkdir=yes`. Delete `outdir`/`dir_exec`. Call `exec_refine3D_workflow(cline_s, l_local)`. The
      child makes its own numbered directory inside the sweep's directory, as a direct run does.
   5. After each child, return to the sweep's directory with `simple_chdir` and reset `CWD_GLOB`. Record the
      child's directory and work project.
   6. After the last child, merge (decision 7), register maps (decision 8), rebind references (decision 9) and
      write the parent once with `write(..., tempfile=.true.)`.
7. **Merge.** This is a new `sp_project` routine. For state `s`, compute the parent rows with
   `select_refine3D_state_rows` on the sweep's unchanged copy. That gives the same rows the child selected,
   because the child prepared its work project from the same bytes. Verify work row `i` against parent row
   `parent_rows(i)` with the stack/index check of `prepare_refine3D_state_work_project` (lines 470-481), extracted
   into an `sp_project` routine that both use (Phase 1).
   - Rows whose parent label is `s` take the child's `ptcl3D` row. In hard mode that is every selected row; in flex
     mode, the hard-labelled rows.
   - The parent's identity keys are then restored: `state`, `stkind`, and every other key Phase 0 lists as identity.
   - The child's `stkind` is never kept, because pruning renumbers stacks.
   - A row whose half or observation fields differ from the parent's is an error, not a relabelling.
   - Flex rows labelled with another state keep the parent's values.
8. **Out segment.** Each child's state-1 `vol` and `fsc` entries, with their `pop`/`mass`/`ess` fields, replace
   state `s`'s entries in the merged parent. Use `remove_state_artifacts_from_osout` and the existing
   `add_vol2os_out`/`add_fsc2os_out`, which store absolute paths. The single `vol_msk` entry keeps the last state's
   mask. The weight-set registration is the parent's, unchanged.
9. **References rebound.** Each child's reference record is valid against its work project, not against the merged
   parent, and `preflight_cont_startup` checks the whole project's layout and the producer state (lines 266-281).
   A new routine in `simple_refine3D_provenance.f90` rebinds a child's record:
   - it reads and validates the record;
   - it keeps the record's constituents, settings digest, map geometry and `align_lp`;
   - it recomputes particle identity on the merged parent and sets the producer state to `s`;
   - it copies the child's recorded state-01 masks (automask, NU envelope, order) to their state-`s` names beside
     the map, and names those copies in the record. There is no copy for `s=1`.

   The record format is unchanged. A child with no record keeps the legacy route. A pending or invalid record fails
   the sweep. A rebound child work project stops being a continuation source; continuation is from the merged
   parent.

   Local staging already copies recorded automasks and NU envelopes from producer state `s` to state 01.
   Extend the same `stage_reference_derivatives` routine in `simple_refine3D_strategy.f90` to do this for
   `ORDER_ENVMASK_FBODY` too, following the existing cleanup of a stale destination when the mask is unrecorded.
10. **Stores and trailing chains.** The sweep does not copy or register the children's sigma2 stores or trailing
    chains in the parent. A later local run validates the parent's sigma2 store as today. Its existing chain
    validator checks any chain beside the source map against the selected child work project, not the full merged
    parent. Reuse a compatible chain; when fractional trailing is active and the chain is absent or incompatible,
    use the existing full-first-update fallback. A hard local→local sweep can recreate a compatible work project,
    so rejection is not guaranteed. No sweep-specific chain reset or change to the validator is needed.
11. **Halves.** The sweep refuses a project in which a row it would refine has no half, and names the reason.
    Children must not each label their own subset.
12. **Interruption.** The merged parent appears only through the final rename. An interrupted sweep leaves the
    sweep's unchanged parent copy and child directories. Rerun from the original input; there is no resume.

Defaults taken without a question:

- a log line per state: state, rows, applied mass and effective sample size under flex, run directory, seconds;
- a table `refine3D_state_sweep.txt` in the sweep's directory (state, child directory, work project), named in
  `simple_refine3D_fnames.f90`;
- ascending state order.

## 4. Not in this plan

- A joint multi-state engine, per-state sigma2 spectra, concurrency between states, and restarting an interrupted
  sweep.
- A per-state pose sidecar for flex rows labelled with another state. Today's direct `state=X` runs seed those rows
  from the parent's poses, and the sweep does no worse. The report gives the measured divergence (section 6); a
  sidecar is a later plan only if real data show it matters.
- Any change to `update_job_descriptions_in_project`, `refine3D_states`, `refine3D_states_cavgs`, `flex_pca` or
  `solve3D`, beyond decision 2's message.
- The ruling of 2026-11-27 that `refine3D_local refine=prob_neigh` requires `cont_view` (unset becomes `lm`, `no`
  refused, the 20,000/500 grid defaults go). It waits for the PorF read of `polar_continuous_view_report.md`,
  section 8, and lands with that read, not here.
- Folders out of bounds: `src/main/strategies/search/**`, `src/main/pftc/**`, `src/main/cftc/**`,
  `src/main/sigma2/**`.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 1 | identical | fast gate; `unit_project`; direct `state=X` identity under the Phase 0 rule |
| 2 | new path, off-state identical | Phase 1's identity; `unit_project` merge test; sweep cases in `flex_pca_blobs` |
| 3 | new path (rebinding and order-mask staging) | provenance and existing strategy testers; local-after-auto sweep cases in `flex_pca_blobs` |
| 4 | docs, report | fast gate |

## 6. Acceptance

- **Fixture.** `flex_pca_blobs` (`simple_commanders_test_highlevel_flex.f90`): two blob states, 2,000 particles, a
  frozen weight set. `run_state_refinements` refines each state directly with `m_estimator=flex` and `no`
  (`STATE_REFINE_MAXITS = 5`). The gate is `STATE_REFINE_MARGIN = 0.01`, measured on the profile of 2026-10-19.
- **Identity rule** (measured, Phase 0). Two direct `refine3D_auto state=X` runs per mode at `nthr=1`, fixed seed.
  If the final maps and the work `ptcl3D` are byte-identical, identity gates are byte gates. Otherwise the tolerance
  is twice the largest difference between the two runs, frozen in Phase 0. Exclusions: paths, `projinfo`,
  timestamps, logs.
- **Merge** (`unit_project`). The fixture has several stacks, drops a whole earlier stack from one state's
  selection, and keeps sparse indices in the surviving stack. Checks:
  - copied fields equal the child's;
  - identity keys, labels and halves equal the parent's;
  - rows with other labels are unchanged;
  - every merged row resolves to its original stack file and index;
  - a foreign or reordered work project is refused.

  The test must fail if the child's `stkind` is kept.
- **Sweep** (`flex_pca_blobs`, Phase 2). Run `refine3D_auto` without `state=` on the published project, once with
  `m_estimator=no` and once with `flex`. Checks:
  - both states run, each in its own directory;
  - each child's map meets the existing truth gate against the direct run of the same state and mode (corr ≥ direct
    − `STATE_REFINE_MARGIN`);
  - the merged parent's labels equal the input's;
  - state `s`'s `vol`/`fsc` entries name child `s`'s files;
  - under flex, the weight set still validates;
  - a relative input-file argument such as `pcg_mskfile` resolves to the entry directory's file in both children
    (include it in one existing sweep case; the matching direct run uses the same file);
  - the command-line project is byte-identical, as the existing check does.

  Once, outside the entry, at `nthr=1`: child `s`'s map and work `ptcl3D` equal a direct run's under the identity
  rule.
- **Chain** (`flex_pca_blobs`, Phase 3). Run `refine3D_local` without `state=` on each merged parent. Checks:
  - it starts without a provenance refusal;
  - for `s ≥ 2` it stages the translated masks, checked by checksum against the child's state-01 originals;
  - each state's local map correlates with truth at least as well as the auto sweep's map, minus
    `STATE_REFINE_MARGIN`.

  The margin was measured for direct auto runs; its use here is guessed. So Phase 0 measures direct auto→local per
  state and mode twice on the base, and stops if the base itself fails this rule.
- **Rebinding** (provenance tester). A record rebound to state 2 validates against the merged project, names the
  state-02 masks, and keeps the child's settings digest and `align_lp`. A pending record is refused.
- **Order-mask staging** (existing strategy tester). A recorded state-02 order mask is staged under state 01
  with identical bytes; an unrecorded mask is not copied and a stale state-01 destination is removed.
- **Off-state.** The existing direct `state=X` runs of `flex_pca_blobs` and the single-state `refine3D_local` runs of
  `cont_refine3D_1jxy` pass unchanged.
- **Reported, not gated:**
  - each state's FSC (`fsc_state01.bin` of each child) and the sweep table;
  - under flex, the median angular and shift difference between the two states' refined poses on the rows that
    both selected. This is the evidence for or against a later pose sidecar.
- **Cost.** Phase 0 records `flex_pca_blobs` wall time. The new cases add four auto children and four local
  children. Stop if the entry exceeds three times its Phase 0 wall time, or 7,200 s; a separate entry is then
  Hans's call.

## 7. Common rules and stop conditions

- No new command-line parameter, no environment knob. No floor, tolerance or assertion relaxed. Fix as you go, with
  a file-table row first. The fast gate runs at one thread; the high-level entry at its registered count, except
  the `nthr=1` identity runs.
- **Stop** if:
  - a child changes a prepared single-state label;
  - the merge would change a parent label, stack, index or half;
  - Phase 0's identity runs disagree by more than `STATE_REFINE_MARGIN`;
  - rebinding needs a change to the record format;
  - an edit is needed outside the file table or in an out-of-bounds folder;
  - the cost rule above is exceeded.
- **No ruling waives:** the fast gate, Phase 1 identity, the merge test, and the sweep and chain cases.

## 8. Phases

- **Phase 0, baseline.** Build Release and Debug, run the fast gate, record the machine and its load. Then:
  - run `flex_pca_blobs` twice, recording wall time and gate values;
  - run the identity runs and the direct auto→local chain reference of section 6;
  - record what `refine3D_auto` without `state=` does today on the fixture's two-state hard project (auto has no
    refusal for it), so the behaviour change is deliberate;
  - list the `ptcl3D` keys a child writes, and freeze the identity keys of decision 7;
  - check the fixture's selected rows all carry halves.

  Exit: Progress holds the identity verdict or tolerance, the chain reference per state and mode, the wall time,
  the base behaviour, and the key lists.
- **Phase 1, row map.** The parent-to-work row check becomes one `sp_project` routine, and
  `prepare_refine3D_state_work_project` calls it. Exit: a `unit_project` test (renumbered stacks and sparse indices
  accepted; a foreign or reordered stack refused); direct `state=X` identity.
- **Phase 2, sweep.** Decisions 4-8 and 10-12:
  - the trigger and the dispatch;
  - `run_refine3D_state_sweep`;
  - the merge and out-entry routines;
  - decision 2's message.

  Exit: the merge test and the sweep cases of section 6; the off-state cases pass.
- **Phase 3, references.** Decision 9, the order-mask extension to existing local staging, and the chain cases.
  Keep decision 10's existing compatible-chain reuse/fallback behavior. Exit: the rebinding and order-mask staging
  tests and the chain cases of section 6.
- **Phase 4, docs and report.** Update the living descriptions: `refine3D_policy.md` ("Per-state M-estimation", a new
  "State sweep"), `refine3D_local_policy.md`, `refine3D_auto_policy.md`, `heterogeneity/refine3D_states_policy.md`
  (decision 2) and `how2_process_heterogeneous_datasets.md`. Write the report `refine3D_state_sweep_report.md`.
  Exit: fast gate; every listed policy names the sweep; the report exists.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | dispatch in the auto/local commanders, input-path resolution before the child template, `run_refine3D_state_sweep`, row-map call, `refine3D_states` message | 1-3 |
| `src/main/project/simple_sp_project.f90`, `simple_sp_project_ptcl.f90`, `simple_sp_project_out.f90` (if needed), `simple_sp_project_tester.f90` | row map, merge, out-entry relabel; tests | 1-2 |
| `src/main/volume/simple_refine3D_provenance.f90`, `simple_refine3D_provenance_tester.f90` | record rebinding with mask translation | 3 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `simple_refine3D_strategy_tester.f90` | extend existing local staging to recorded order masks; recorded/unrecorded mask tests | 3 |
| `src/defs/simple_refine3D_fnames.f90` | sweep table filename | 2 |
| `src/main/ui/simple/simple_ui_refine3D.f90` | help of `state` and `m_estimator` | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_flex.f90` | sweep and chain cases | 2-3 |
| `doc/policies/3D/refine3D_policy.md`, `refine3D_local_policy.md`, `refine3D_auto_policy.md`, `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/how2s/how2_process_heterogeneous_datasets.md` | living descriptions | 4 |
| `doc/implementation_notes/planned/refine3D_state_sweep*.md` | this plan, Progress, report | all |

## Progress

2026-11-27: proposal only; no implementation phase has started. Rewritten after ruling 3 from a 434-line revision
that had grown through two reviews. That revision's review is `refine3D_state_sweep_review.md`.

The rewrite drops:

- the state pose sidecar with its continuation payload;
- the change to the CLI job-description write;
- the fault-injection interruption matrix;
- the hard/flex × gridding/PCG fractional matrix;
- the in-memory staging across children;
- the `cont_view` default switch, moved to section 4.

Still covered from the earlier reviews:

- child stack renumbering (decision 7, merge test);
- re-entry with a normalized command (decision 5: dispatch before mutation);
- earlier merges changing later seeds (moot: one merge at the end);
- state-01 masks under relabelling (decision 9);
- partial publication (one atomic write);
- prepared-versus-parent labels (Phase 0).

2026-11-27 follow-up: resolve relative input paths before saving the child template; include order-mask staging
in the existing strategy routine and file table; leave compatible-chain reuse and the missing/incompatible-chain
fallback to the existing validator on the selected work project. The simplified design is unchanged.
Documentation only; no implementation, build or runtime tests in this revision.
