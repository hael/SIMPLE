# refine3D_states_cavgs: design note

**Date:** October 8, 2026
**Status:** Proposal for review, revised the same day after an independent
review (Codex). Written from a read of GitHub master at commit `89bbc49`.
Nothing was compiled or run.

## Summary

`refine3D_states_cavgs` is a proposed program that sorts the sub-class averages
written by `cls_expansion` into conformational states. It runs the
`refine3D_states` frequency march on a temporary class-average work project,
the way `solve3D_cavgs` runs `solve3D`.

The program uses no nonuniform (NU) machinery and no automasking. It uses the
FSC between the even and odd class-average maps the way `solve3D_cavgs` does.
Its product is a state label and a pose per sub-class average in `cls3D`,
handed back to the particles, plus class-average state maps.

Its natural role is a fast state initializer beside `flex_pca`. The
particle-level `refine3D_states` then continues from populated labels and owns
the final particle maps.

## Scope and non-goals

The change adds one public program, a controller extracted from the states
commander, a shared work-project builder, one option on the pose initializer
and two project routines; the base `refine3D` code is not touched.

In scope:

- A new program `refine3D_states_cavgs`, wired through ui, exec and commander.
- Input: a project after `cls_expansion` that also holds a consensus map
  (`vol`, state 1) and particle poses from a prior 3D refinement.
- Output: state, pose and score per sub-class average in `cls3D`; a state label
  per particle in `ptcl3D`; class-average state maps registered as `vol_cavg`.

Not in scope:

- Running the base `refine3D` on the `cls3D` segment. The march runs on
  `ptcl3D` rows of a temporary work project.
- Final particle maps. Those stay with particle-level `refine3D_states`; this
  program ships class-average maps, as `solve3D_cavgs` does.
- `flex_pca` initialization, the weighted mode (`m_estimator`), NU filtering
  and automasking.
- Any change to `cls_expansion`, or any behaviour change in `solve3D_cavgs` and
  particle `refine3D_states`.
- Moving `classify3D_refs` onto the shared controller. It can follow later.

## Work project, as in `solve3D_cavgs`

Decided on 2026-10-08 with Afan: `refine3D_states_cavgs` analyses a temporary
work project set up exactly as in `solve3D_cavgs`. `cls_expansion` will
generate even/odd class-average pairs in the same form as the class averages
that `solve3D_cavgs` takes as input.

The work project follows section 3 of `solve3D_cavgs_policy.md`:

- It preserves project info, compute environment and job-process metadata.
- It registers the even and odd class-average stacks as two stacks.
- It expands each class average into an even and an odd `ptcl3D` row without
  CTF, carrying class and state.
- It owns a transient canonical sigma2 state file and never inherits the input
  project's.
- It is deleted, with that file, when the run completes.

Everything downstream then sees an ordinary `ptcl3D` lineage, so Euclidean
matching, `reconstruct3D` and `calc_final_rec` run unchanged.

Only the construction is shared with `solve3D_cavgs`. The commit is not:
`solve3D_cavgs` recreates `cls3D`, which here would erase the `cluster` and
`accept` fields that `cls_expansion` wrote.

## Evidence contract: FSC as in solve3D_cavgs, no NU

Decided on 2026-10-08: the program supports FSC the way `solve3D_cavgs` does
and keeps all NU machinery off. With even/odd class-average pairs from
`cls_expansion`, each state gets one half map from the even member averages and
one from the odd.

| What the even/odd FSC may drive | Where | Setting |
| --- | --- | --- |
| ML regularization of the state maps | `add_invtausq2rho(fsc)` in `restore_gridding_pair`, `simple_commanders_rec_distr.f90` | `ml_reg=yes`, from the first block |
| Final maps through the ML route, with canonical sigmas reused or bootstrapped | `calc_final_rec` in `simple_final_rec.f90` | follows `ml_reg` |
| Low-pass of the final snapshot at FSC 0.143 | `solve3D_state_fsc_lowpass`, called by `write_final_rec_outputs` | never finer than `lpstop` |
| FSC and resolution records per iteration | volume assembly | always, diagnostic |

| What stays off | Where it would enter | Closed by |
| --- | --- | --- |
| NU filter bank, NU masks and low-pass promotion | `filt_mode=nonuniform*`: `nonuniform_filter_state`, `refresh_matching_lp_from_project` | `filt_mode=none`, fixed |
| Matching low-pass taken from the FSC | `filt_mode=fsc` or `uniform` sets `l_lpauto`, `simple_parameters_phases.f90` | `filt_mode=none`; the march sets `lp` per block |
| Automasking and envelope-masked FSC | `automsk`, `envfsc` | both `no`, fixed |
| Postprocessing of the final maps | `postprocess_states` in `simple_rec3D_service.f90` | `l_postprocess=.false.` |

The second table is what `solve3D_cavgs` already fixes for its class-average
route: `filt_mode=none`, `automsk=no`, `envfsc=no` and no postprocessing.

With `ml_reg=yes` the matching reference is the merged state map, regularized
at assembly, and the block's `lp` limits the matching band.
`simple_refine3D_stage_plan` needs only box, sampling, `lpstart` and `lpstop`.

The FSC is not gold-standard, as in `solve3D_cavgs`: both halves are matched
against the merged reference. With class expansion as it stands today, the
halves also share the fitted basis, the labels, the in-plane alignment and the
Wiener prior, so the curve will read optimistic.

Avoiding a resolution claim does not remove that bias. The FSC is converted to
an SSNR that enters the reconstruction denominator, so an optimistic curve
under-regularizes the state maps.

Decided on 2026-10-08: ML regularization is on from the first block all the
same, because the march starts from good orientations. The comparison against
`ml_reg=no` stays in the validation list as a check on that choice.

Sibling sub-averages do not share particles, because the delivered stack is
restored from hard labels.

## Workflow

The program builds a class-average work project, seeds poses and labels in it,
runs the shared states controller there, and writes labels, poses and
class-average maps back in one final step.

```text
INPUT PROJECT                                  TEMPORARY WORK PROJECT

cls_expansion output            -- build -->   2N rows, even and odd averages
  sub-average stack, even/odd                    no CTF, own transient sigma2 state
  cls2D selection, ptcl2D class                          |
  consensus vol, ptcl3D poses                  seed poses and labels
                                                 one cc pass against the consensus map
                                                 paired, balanced random labels
  (not modified while                                    |
   the march runs)                             refine3D_states frequency march
                                                 lp set per block, even/odd FSC, no NU
                                                         |
                                               gates: half agreement, state support
                                                         |
written back                    <-- commit --  state maps at native sampling
  cls3D: state, pose, score                      no postprocessing, low-pass snapshot
  ptcl3D: state only, poses kept
  out: vol_cavg state maps
```

The input project is read once at build and written once at commit; every FSC
the march computes stays in the work project.

1. **Preflight.** The run stops before any write unless all of these hold:
   - the sub-average stack and its `_even` and `_odd` companions exist and hold
     the same number of images as `cls2D` and `cls3D` have rows;
   - the three stacks and the consensus `vol` of state 1 share box and
     sampling;
   - `ptcl3D` is congruent with `ptcl2D` and carries usable poses from the
     prior refinement, and the two active masks agree;
   - `nstates >= 2`, and the selected sub-averages are enough to support every
     state in each half.
2. **Snapshot the particle active mask.** The `ptcl3D` states of the input
   project are recorded before anything else, for the hand-back in step 10.
3. **Build the work project** with the shared builder. The even and odd stacks
   become two stacks and `2N` `ptcl3D` rows without CTF. Row selection follows
   the `cls2D` state flags. The project gets its own transient sigma2 state
   file. All children run from a separate work command line with the work
   project file and `mkdir=no`; the public command line is never mutated.
4. **Seed poses.** One fixed-reference `cc` pass aligns every row against the
   consensus map at `lpstart`, through
   `initialize_poses_against_external_references` in its new strict-full mode.
   The same service bootstraps sigma2 from the images. After its validation
   the `sampled` and `updatecnt` fields are cleared and written, so the
   coverage check in step 6 counts the march and not the seed pass.
5. **Seed labels in pairs.** Balanced uniform random labels are drawn over the
   `N` sub-averages and copied to both rows of each pair, so every state is
   balanced in both halves. Multi-state labels already present in `cls3D` are
   copied to both rows instead.
6. **Run the shared states controller** with the class-average policy. It
   reconstructs the startup state maps from the labels, runs the frequency
   march under `pose_policy` in strict full-batch mode, and requires every
   active row to be updated by the march.
7. **Apply the gates.** Pairs whose even and odd rows end in different states
   are set to state 0 in the work project. Every requested state must then
   keep at least `min_state_frac` of the committed sub-averages. If a state
   falls short the run aborts, and the input project is untouched.
8. **Reconstruct the final maps** from the rows that passed the gates, with the
   ending `solve3D_cavgs` uses: `calc_final_rec` at native sampling without
   postprocessing, then the raw and low-pass snapshot maps.
9. **Commit to `cls3D` in place.** State, pose and score of each committed
   sub-average are written into the existing rows. The `cluster` and `accept`
   fields from `cls_expansion` are kept. The even/odd state agreement and
   angular distance are logged.
10. **Hand back and clean up, as one final write.** Particles that were
    inactive stay at state 0. Active particles of a committed sub-class take
    its state; `ptcl3D` poses are left as they are. The `vol_cavg` entries are
    replaced as a set, and the work project and its sigma2 state are deleted.

Particle-level `refine3D_states` can then be run on the same project. It finds
populated multi-state labels, reconstructs state maps from them and proceeds
without `flex_pca`.

| Setting | Value | Status |
| --- | --- | --- |
| `filt_mode` | `none` | fixed, anything else is refused |
| `ml_reg` | `yes` | from the first block; `ml_reg=no` turns it off |
| `envfsc`, `automsk`, `combine_eo` | `no` | fixed |
| `objfun`, `sigma_est` | `euclid`, `global` | fixed, as in `solve3D_cavgs` |
| `rec_backend` | `gridding` | fixed, as in `solve3D_cavgs` |
| `nsample` | all active rows | fixed, user values are refused |
| `update_frac`, `cohort_sampling` | unset | fixed, user values are refused |
| `trail_rec` | `no` | fixed |
| `pose_policy` | `global` | overridable, `local` or `global` |
| `nstates` | none | required, at least 2 |
| `min_state_frac` | 0.1 | overridable, used by the support gate |
| `lpstart`, `lpstop` | to be decided | overridable |
| GUI publication | from the input project only | fixed |

All of these live in the typed controller policy, not in command-line state.

## Refactoring steps

Six code changes carry the program. The first two are extractions that must
leave particle `refine3D_states` and `solve3D_cavgs` unchanged, and they land
first.

1. **Shared states controller.** The march now lives in contained procedures
   of `exec_refine3D_states`. It moves into a module, working name
   `simple_refine3D_states_controller`, with a typed policy:
   - the controller owns sampling and frequency planning, state-count
     resolution from labels, startup state maps, stage execution, the coverage
     read, the missing-update pass and the state-overlap read;
   - the policy carries every fixed setting, the update scheme (sampled or
     strict full batch) and the GUI switch;
   - each commander keeps what is specific to its lineage. The particle
     commander keeps `flex_pca` initialization and the postprocessed ending.
     The class-average commander keeps build, seeding, gates, finalization and
     commit.

   Large objects such as `parameters` and `cmdline` are allocatable components,
   as the compile-time policy requires.
2. **Shared work-project builder.** A new module, working name
   `simple_cavgs_work_project`, takes the construction out of
   `exec_solve3D_cavgs`: stack lookup, even/odd companions, project info copy,
   transient sigma2 state path and row expansion (roughly lines 164 to 227).
   The diagnostics `conv_eo` and `conv_eo_states` and the `nparts` cap of
   `configure_cavgs_distributed_clines` move with it. `solve3D_cavgs` keeps its
   own commit and its `map2ptcls`.
3. **Strict-full option on the pose initializer.** An optional argument of
   `initialize_poses_against_external_references` sets the sample to all
   active rows, lifts the 100,000-row cap and requires every row to be
   updated. Without it the three existing callers behave as before.
4. **New commander `commander_refine3D_states_cavgs`.** It owns workflow steps
   1 to 5 and 7 to 10.
5. **Two project routines** in `simple_sp_project_ptcl.f90`: an in-place
   update of state, pose and score in existing `cls3D` rows, and a state-only
   hand-back that respects a supplied active mask and leaves poses alone.
6. **Registration** in the UI and the exec router.

| File | Change |
| --- | --- |
| `simple_refine3D_states_controller.f90` (new, location to confirm) | controller and typed policy |
| `simple_commanders_refine3D.f90` | `exec_refine3D_states` calls the controller; new commander |
| `simple_cavgs_work_project.f90` (new, location to confirm) | work-project builder, even/odd diagnostics |
| `simple_commanders_solve3D.f90` | `exec_solve3D_cavgs` calls the builder |
| `simple_external_reference_pose_initialization.f90` | optional strict-full argument |
| `simple_sp_project.f90`, `simple_sp_project_ptcl.f90` | in-place `cls3D` update; state-only hand-back |
| `simple_ui_heterogeneity.f90` | `new_refine3D_states_cavgs` |
| `simple_exec_refine3D.f90` | router case `refine3D_states_cavgs` |

Alternative rejected: entering `exec_refine3D_states` through an in-process
command-line key, as `solve3D_addon` enters `exec_solve3D`. It needs less code
movement, and the first draft of this note leaned that way. The review points
in this revision change the balance. The class-average path now differs in
update scheme, counter handling, gates, ending and commit, which is too much
to hide behind a mode branch and a key outside the argument vocabulary.

The cost is real: about 600 lines of contained procedures over host variables
move, and the particle path must stay identical. `exec_classify3D_refs`
carries second copies of the coverage and missing-update helpers and can adopt
the controller later.

## Open decisions

Six choices remain open; a lean is given for each.

| Decision | Options | Lean |
| --- | --- | --- |
| Half disagreement at the end | Agreement required; a joint decision from both rows' scores; or best half wins, as `solve3D_cavgs` | Agreement required. Disagreeing pairs leave before the final reconstruction, so `vol_cavg` describes exactly the handed-back sub-classes. |
| Snapshot low-pass | FSC 0.143 as in `solve3D_cavgs`; or the same, never finer than `lpstop` | Never finer than `lpstop`, because the FSC reads optimistic. |
| Failed support gate | Abort; or restart with new labels, as the conditional restarts of `solve3D_cavgs` | Abort in the first version, restarts later. |
| Active particles without a committed sub-class | Excluded at state 0; or classified later against the state maps | Excluded, with the count logged and the original mask kept so a later pass can re-admit them. |
| Label balance | Pairs balanced globally; or also across orientation bins | Global first. Add orientation bins if the phantom shows a view bias. |
| Frequency band | `refine3D_states` uses 10 to 6 A; `solve3D_cavgs` uses 20 to 8 A, or limits derived from class FRCs | 20 to 8 A to start. The covariance fit behind the sub-classes stops at 30 A by default. FRC-derived limits need sub-class FRCs in the project. |

Settled in this revision: even/odd rows, the shared controller, pair seeding,
strict full-batch operation, the counter reset after pose seeding, the in-place
commit, inactive particles staying inactive, and ML regularization from the
first block.

Active particles without a committed sub-class come from three sources: parents
too small for `cls_expansion` to split, deselected sub-classes, and pairs that
disagreed. Class 0 in `ptcl2D` also holds particles that were already inactive,
which is why the active mask is snapshotted first.

Two smaller points:

- Reference low-pass shape. `refine3D_states` uses the cosine edge at `lp`;
  `solve3D_cavgs` uses `gauref=yes`. This applies only when ML regularization
  is turned off. The cosine edge keeps the two states programs matched.
- Per-particle labels are a seed, not a result. The `cls_expansion` note states
  that the latent orders members but does not label individual particles. The
  particle-level march must stay free to move them, which it is.

Current behaviour of `cls_expansion` that bears on this program, confirmed by
two independent reads of the code. The final restoration pass runs with labels
taken from the project, so the weights are one-hot and the separation is zero.
It rewrites `cls_expansion_weights.txt` and rebuilds `cls2D`. The delivered
`sep` is therefore 0, `neff` equals `pop`, and the soft responsibilities are
gone. That matters if `sep` should gate which sub-averages enter, or if a soft
hand-back is wanted later.

## Validation and documentation

Compilation and runtime tests are performed by the user; nothing here has been
built.

The primary automated test builds on the generated two-state FLEX phantom
(`create_flex_pca_phantom_fixture` and `build_flex_pca_phantom_project` in
`simple_flex_pca_application_tester.f90`). The EMD-8440/8445 mix used for the
`cls_expansion` benchmarks stays as a secondary, manual benchmark.

Regression, before the new program exists:

- [ ] Particle `refine3D_states` gives the same outputs before and after the
  controller extraction.
- [ ] `solve3D_cavgs` gives the same outputs before and after the builder
  extraction, single-state and multi-state.
- [ ] The three existing callers of the pose initializer are unchanged without
  the strict-full argument.

Assertions on the phantom:

- [ ] Full-batch coverage: every active row is updated in every iteration, and
  no fractional or trailing flag is set.
- [ ] Counter reset: `sampled` and `updatecnt` are empty between pose seeding
  and the first march iteration.
- [ ] Metadata preservation: `cluster` and `accept` in `cls3D`, and all of
  `cls2D`, are unchanged by the commit.
- [ ] Half disagreement: a forced disagreement removes the pair from the final
  maps and from the hand-back.
- [ ] Collapse: a forced empty state aborts the run and leaves the input
  project unchanged.
- [ ] Inactive rows: particles at state 0 before the run are at state 0 after
  it.
- [ ] Shared-memory and distributed runs give equivalent labels, including
  `nparts` larger than the row count.
- [ ] State recovery: purity of the committed sub-averages and of the
  handed-back particle labels against the planted truth, with and without ML
  regularization.

Further checks:

- [ ] The policy refuses `filt_mode=nonuniform_lpset`, `automsk=yes`,
  `update_frac`, `nsample` and `trail_rec=yes` before the first iteration.
- [ ] No NU product appears: the run directory holds no `_nu_filt` volumes and
  no NU envelope masks.
- [ ] End to end: particle `refine3D_states` on the output project starts from
  the handed-back labels without `flex_pca`.
- [ ] The input project's out segment keeps its `vol` and `fsc` entries, and
  its `vol_cavg` entries are replaced as a set or not at all.
- [ ] A failure at any step leaves the input project unchanged and removes the
  work project and its sigma2 state.

Documentation to add or update once the design is agreed:

- New `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md`, laid out
  like `solve3D_cavgs_policy.md`.
- `refine3D_states_policy.md`: name `refine3D_states_cavgs` as a second source
  of populated labels beside `flex_pca`, and the controller as the owner of the
  march.
- New page in `doc/algorithms/heterogeneity_analysis/` and a line in its
  README.
- `doc/how2s/how2_process_heterogeneous_datasets.md`: the route
  `cls_expansion`, then `refine3D_states_cavgs`, then `refine3D_states`.
- The `simple-refine3d` skill code map, only with approval.
