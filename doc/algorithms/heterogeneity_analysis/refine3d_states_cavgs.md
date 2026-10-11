# Class-Average State Refinement (`refine3D_states_cavgs`)

## Problem

Given the sub-class averages that `cls_expansion` splits from each 2D class,
and a consensus map, sort the sub-averages into conformational states and
hand the states back to the particles. Class averages carry far more signal
per image than particles, so the states can be found and refined on a few
hundred images before particle-level
[`refine3D_states`](refine3d_states.md) continues from the handed-back labels.
The shared estimator and bandwidth schedule are described in the
[heterogeneity-analysis overview](README.md).

## Algorithm

Each sub-average enters a temporary work project as two independent rows, its
even and its odd average, without CTF, as in `solve3D_cavgs`. For `N`
sub-averages the rows are:

```text
rows 1..N      even averages   (eo = 0)
rows N+1..2N   odd averages    (eo = 1)
```

The run then proceeds in five steps:

1. One fixed-reference cross-correlation pass aligns every row against the
   consensus map at `lpstart` and reconstructs a checkpoint map from the class
   averages themselves.
2. [`flex_pca`](flex_pca.md) maps the rows to at most `nstates` discrete
   states at fixed poses, with the checkpoint map as its mean, its covariance band at
   `lpstart` and the native box. Labelled input (two or more populated `cls3D`
   labels) skips this step.
3. The `refine3D_states` frequency march refines poses and states, with every
   row updated in every iteration and one common band from `lpstart` to
   `lpstop`.
4. For each sub-average, the row with the better correlation gives the pose,
   state and score. Disagreement between the halves is measured, as the
   even/odd state overlap and angular distance, but not resolved.
5. State maps are reconstructed from all rows at native sampling. The
   populated states are numbered consecutively and written back: per
   sub-average into `cls3D`, per particle (state only) through each particle's
   2D class.

Every reconstruction uses the PCG backend by default, on the spherical
`mskdiam` support.

## Rationale

Class averages are sparse but clean data, so a regularized iterative
reconstruction (PCG) suits them. Treating the even and odd averages as
independent rows keeps the half maps separate and copies what
`solve3D_cavgs` does successfully. The FSC still reads optimistic: both halves
are matched against the merged reference, and the halves share the class
expansion that made them. The low-pass snapshot is therefore never finer than
`lpstop`.

Flex initialization replaces random state seeding. The march then lets poses
and states move together, as in the particle workflow. An empty state is
reported and dropped rather than failing the run, and the consecutive
renumbering lets particle `refine3D_states` start from the labels directly.

## Implementation

- Workflow: `exec_refine3D_states_cavgs` in
  `src/main/commanders/simple/simple_commanders_refine3D.f90`; flex
  initialization and the march in `src/main/simple_refine3D_states_controller.f90`.
- Work project and even/odd diagnostics: `src/main/simple_cavgs_work_project.f90`.
- Pose seeding: `src/main/simple_external_reference_pose_initialization.f90`.
- Hand-back: `map_cls3D_states2ptcls` in `src/main/project/simple_sp_project_ptcl.f90`.
- Policy: `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md`.
