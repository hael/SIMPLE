# Local 3D Refinement Policy

The command boundary and startup/fractional source paths are implemented.
Compilation and runtime validation for the current changes are pending on Linux.

`simple_exec prg=refine3D_local projfile=<aligned.simple> mskdiam=<A>`
continues single-state Euclidean refinement from stored particle poses using
`refine=cont`. The command replaces `refine3D_auto pose_cont=only`. Auto retains
`pose_cont=no|yes`; base `refine3D refine=cont` remains available.

Before parameter initialization, the local command requires an existing merged
map and its even/odd pair, compatible dimensions/sampling/physical extent,
valid saved poses and half labels, and selected particles in both half-sets.
Producer records bind the reference constituents to their project/state,
geometry, generation, particle identity and settings. Missing records admit a
legacy base pair with a logged limitation and FSC regeneration. Incomplete or
inconsistent recorded sets are rejected before parameter initialization.
Recorded continuation inputs also require compatible ordered rows, half labels,
active population and selected producer state. Explicit external references from
a different particle lineage retain the caller-supplied alignment/half contract;
same-lineage half or row changes cannot use that exception. These checks do not
establish independent reconstruction or a common coordinate frame.

Explicit `vol1` takes precedence over the project's registered volume. A
multistate project requires `state=X`. The local and auto commands share state
work-project preparation, including frozen `m_estimator=flex` row selection.
The local preflight checks the same selected rows before parameter setup.

The local command rejects `refine`, `pose_cont`, projection-grid arguments,
discrete pose acquisition and registration, and overrides that change Euclidean
scoring, ML regularization, fixed reference centering or separate halves.
Local rejects explicit `lp` and `nonuniform_lpset`; use `lpstop` to cap the band
without merging reference halves. Local controls include `athres_cont`, `trs`, `cont_conv_rot`, `cont_conv_shift`,
`cont_conv_frac`, `nsample`, iteration limits and the existing backend controls.
An explicit `trs` is retained; an absent value uses the existing derived bound.

The commands share one workflow body and one UI input builder. Local captures
source paths before run-directory creation, stages the saved map pair, reuses a
committed sigma file only after validation against the selected work-project
rows and its producer record, and delegates a missing or incompatible sigma state to the existing
bootstrap owner. It skips the startup particle-power and reconstruction passes.
A missing FSC is calculated from the saved half maps. The matcher prepares its
Cartesian references from those maps. A compatible recorded NU bundle carries
its filtered maps, local-resolution map, recorded masks and matching LP into the
run. Missing or incompatible derivatives use the existing logged FSC-filtered
fallback; volume assembly writes updated FSC, mask and NU products.

A compatible complete trailing chain is copied as a set with its manifest last.
Without one, local makes iteration 1 a full update and seeds the chain; later
iterations use the configured fraction. Shared and distributed strategies use
the same first-update transition. The shared `calc_final_rec` ending remains in
place. Local defaults to `autoscale=no` (native sampling); auto retains
`autoscale=yes`. Explicit settings take precedence.

The remaining design and Linux acceptance checks are in the living
[local-workflow plan](../../refactoring_notes/planned/refine3d_local_workflow.md).
