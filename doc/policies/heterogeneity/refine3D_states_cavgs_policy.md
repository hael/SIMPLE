# Refine3D States Cavgs Policy

This document records the policy for `refine3D_states_cavgs`, which sorts the
sub-class averages written by `cls_expansion` into conformational states. It is
the class-average counterpart of particle-level `refine3D_states`
([refine3D_states_policy.md](refine3D_states_policy.md)) and builds its work
project as `solve3D_cavgs` does
([solve3D_cavgs_policy.md](../3D/solve3D_cavgs_policy.md)). The design record
is `doc/implementation_notes/completed/refine3D_states_cavgs.md`.

Primary implementation:

- `src/main/ui/simple/simple_ui_heterogeneity.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
  (`exec_refine3D_states_cavgs`)
- `src/main/simple_refine3D_states_controller.f90` (the `flex_pca` state
  initialization and the frequency march, shared with `refine3D_states`)
- `src/main/simple_cavgs_work_project.f90` (the work project and the even/odd
  diagnostics)
- `src/main/simple_external_reference_pose_initialization.f90` (pose seeding)
- `src/main/project/simple_sp_project_ptcl.f90` (`map_cls3D_states2ptcls`)

## 1. Scope

`refine3D_states_cavgs` is a careful analysis of the dynamics that
`cls_expansion` resolves, aimed at good class-average state maps. It is not an
initial-model step. Its products are:

- a state, pose and score per sub-class average in `cls3D`;
- a state label per particle in `ptcl3D`, poses untouched;
- class-average state maps, registered as `vol_cavg`.

Particle-level `refine3D_states` then continues from the handed-back labels
without `flex_pca` and owns the final particle maps.

The program is visible at the Developer level of the user interface. It uses
no nonuniform (NU) filtering and no automasking.

## 2. Inputs

The input is a project after `cls_expansion` that also holds a consensus map
(`vol`, state 1). The run stops before any work unless all of these hold
(`N` is the number of selected sub-averages):

- the class-average stack (registered as `cavg`) exists and holds one image per
  `cls2D` row, `cls3D` has as many rows, and the stack has the box and sampling
  of the particles;
- the `_even` and `_odd` companions exist and match the merged stack in box,
  image count and header sampling;
- the consensus map exists at the box and sampling of the particles;
- `ptcl2D` and `ptcl3D` have the same size, and every particle active in
  `ptcl3D` is active in `ptcl2D`;
- `mskdiam` is given: it is the support of the PCG solves;
- at least 100 sub-averages are selected in `cls2D`;
- for the `flex_pca` initialization (section 6): `nstates >= 3` and flex's
  population floor fits, `nstates * max(1, nint(min_state_frac * R)) <= R`, `R`
  the active rows (`2N` less the halves `cls2D` flags with state 0);
- with populated `cls3D` labels (section 6), an explicit `nstates` equals
  their number.

On the `flex_pca` path, each of flex's two fits must also get at least 100
rows, counted with `flex_pca_half_of`; that check runs right after the work
project is built and before any refinement. Particle poses are not read.

## 3. Fixed Settings and Defaults

The route fixes:

- `filt_mode=none`, `automsk=no`, `envfsc=no`, `combine_eo=no`
- `objfun=euclid`
- `trail_rec=no`, `balance=none`, `greedy_sampling=no`, `frac_best=1`
- `lplim_crit=0.5`, `incrreslim=no`
- `autoscale=yes` (the frequency planner owns downscaling)

A different user value of `filt_mode`, `automsk`, `envfsc`, `combine_eo`,
`trail_rec`, `objfun` or `pcg_solvent` is refused, and so are
`m_estimator` other than `no`, `nsample`, `update_frac`, `cohort_sampling`,
`trail_seed`, `pcg_mskfile`, `vol1..volN`, `multivol_mode` and
`prob_neigh_mode`.

When unset, it supplies:

- `rec_backend=pcg` (`gridding` is accepted as an override and then applies to
  every reconstruction of the run, `flex_pca` included)
- `pose_policy=global`
- `ml_reg=yes`
- `lpstart=10`, `lpstop=6`
- `min_state_frac=0.1`
- `overlap=0.99`
- `prob_inpl=yes`, `center=no`, `keepvol=no`, `mkdir=yes`

`pose_policy` and the local pose bounds behave as in `refine3D_states`.

## 4. Work Project

The workflow writes a temporary project named
`refine3D_states_cavgs_tmpproj.simple`, built as `solve3D_cavgs` builds its own.
The work project:

- preserves project info, compute environment and job-process metadata;
- registers the even and odd class-average stacks as two stacks;
- expands each sub-average, selected or not, into an even row and an odd row
  without CTF: for `M` sub-averages, rows `1..M` are the even averages
  (`eo=0`) and rows `M+1..2M` the odd ones (`eo=1`);
- takes the selection from the `cls2D` states and mirrors it into `ptcl2D`, so
  the rows of a deselected sub-average are inactive;
- makes a row inactive when `cls2D` flags its half with state 0
  (`state_even`, `state_odd`; absent means active), and a sub-average with one
  active row takes that row (section 8). `cls_expansion` writes the flags: 0
  for a sub-class half without members, whose average is only the class's half
  mean, 1 otherwise;
- owns a transient canonical sigma2 state,
  `refine3D_states_cavgs_sigma2_state.bin`, and never inherits the input
  project's;
- copies no out segment.

If distributed execution requested more partitions than rows, `nparts` is
reduced to `2M`. The work project, its sigma2 state and the partition stacks
are deleted when the run completes. A failed run leaves them in the run
directory and the input project unchanged.

## 5. Pose Seeding

One fixed-reference `objfun=cc` pass aligns every active row against the
consensus map at `lpstart`, through
`initialize_poses_against_external_references` with `l_all_rows`: no cohort
cap, every active row must be updated, and `sampled` and `updatecnt` are
cleared after the checkpoint so the coverage check counts only the march. The
pass bootstraps sigma2 from the images and reconstructs a checkpoint map from
the aligned rows. It is iteration 1 of the run.

## 6. State Initialization

Without multi-state labels, `flex_pca` maps the rows to discrete states. It
runs once, on the work project, with the settings of the class-average
policy:

| Setting | Value |
| --- | --- |
| `rec_backend`, `rec_states_backend` | the run's backend, `pcg` by default |
| `box_crop` | the native box |
| `lp` | `lpstart` |
| `vol1` (mean map) | the pose-seeding checkpoint map, copied to `refine3D_states_cavgs_mean.mrc` |
| `npreimages` | `nstates` |
| `min_state_frac` | the run's value, 0.1 by default |
| `min_neff` | 20 |
| distribution | shared memory |

Flex keeps its own cross-validation; its code is not changed. Rows of a state
flex drops go to state 0. Its state weight set is discarded: the march works on
hard labels. The number of states flex delivers becomes the run's `nstates`.

When the selected sub-averages carry two or more populated `cls3D` labels, flex
is skipped. The labels are numbered consecutively (an empty label is dropped),
copied to the active rows of each sub-class (a half flagged with state 0 stays
inactive), and the startup state maps are reconstructed from them. Flex's
preflight counts the active rows the same way, so a flagged half never counts
towards a fit.

## 7. Frequency March

The shared states controller runs the `refine3D_states` frequency march
(`refine=prob_neigh`, one common band from `lpstart` to `lpstop`) with these
differences:

- every active row is updated in every iteration (strict full batch): no
  fractional update, no cohort sampling, no trailing reconstruction;
- the iteration count continues after the pose seeding;
- without an explicit `maxits` the march runs at most 10 iterations (an
  explicit `maxits` must be at least 5);
- no stage maps are published to the GUI.

With the low-pass set per block, both halves are matched against the merged
state map. The FSC compares the half map from the even rows with the one from
the odd rows and is not gold-standard, as in `solve3D_cavgs`.

A missing-update pass then updates any active row the march did not reach.

## 8. Better Row

For each sub-average, the active row with the better correlation gives the
pose, state and score; a tie goes to the odd row, and a sub-average with only
one active row takes that row. A sub-average with no active row gets state 0.
The even/odd state overlap and average angular distance are logged.

The states that hold committed sub-averages are then numbered consecutively,
and the mapping is logged. An empty state never fails the run: it is reported
and dropped. If only one state stays populated, the run completes and says so.

## 9. Final Maps

Final maps are reconstructed by `calc_final_rec` from all active rows as
assigned, at native sampling and without postprocessing. For every populated
state the run writes `rec_final_stateNN.mrc`, `NN` the consecutive state number,
and a low-pass snapshot at the FSC 0.143 point, never finer than `lpstop`
because the even/odd FSC of class averages reads optimistic.

## 10. Commit and Hand-Back

The input project (under `mkdir=yes`, its copy in the run directory) is read
once at the start and written once at the end:

- `cls3D` is updated in place: each committed sub-average's state, Euler
  angles, shift, correlation and projection index are written into its row,
  and the `cluster` and `accept` fields from `cls_expansion` are kept;
- `ptcl3D` states are handed back over the `ptcl2D` selection, snapshotted at
  the start: a selected particle whose sub-average was committed takes its
  state, every other particle gets state 0, and poses are left as they are;
- the `vol_cavg` entries of the out segment are replaced as a set; the other
  entries, the consensus `vol` among them, are kept.

Because the hand-back covers the `ptcl2D` selection, a rerun on the output
starts. It takes the labels path of section 6, so a sub-average the earlier
run left at state 0 stays out, and so do its particles, unless its `cls3D`
state is reset first.

## 11. Validation

User-side validation must cover:

- the `flex_pca` initialization on the class averages (both fits with at
  least 100 rows, the PCG backend, the native box, `lp=lpstart`, the
  checkpoint map as mean) and the labels path;
- full-batch coverage and the counter reset after pose seeding;
- the better-row pick and the logged even/odd overlap;
- preservation of `cluster` and `accept` in `cls3D` and of all of `cls2D`;
- a forced empty state, which must not fail the run;
- particle `refine3D_states` starting on the output from the handed-back labels;
- state recovery against the planted truth on the phantom, with both backends.

The automated gate is a stage of the registered high-level test `flex_pca_blobs`
(`simple_test_exec test=flex_pca_blobs`):
`cls_expansion` and this program on a three-state phantom of 3 600 particles
in 60 exact views, checked against the planted states and directions
(`simple_commanders_test_highlevel_flex.f90`).

No Linux or BOX result is recorded as passing without observed output.
