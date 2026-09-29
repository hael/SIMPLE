# Refine3D Pose-Continuous Workflow Policy

`refine3D_pose_cont` is the automated, single-state workflow that combines
the startup, greedy registration, probabilistic refinement, and final
reconstruction policy of `refine3D_auto` with mandatory Cartesian pose
refinement. It delegates every numerical refinement iteration to base
`refine3D`; it is not a matcher implementation.

## Public modes

The program selects one `pose_cont_mode`. Its public UI defaults to
`post_matcher`. The global parameter default remains `off`, and non-`off`
values are rejected by other programs.

The public `objfun` choices are `euclid` and `cc`. This selects the outer
probabilistic-matcher objective. Cartesian pose refinement continues to use
its internally owned normalized-correlation objective for both choices.

- Startup reconstruction uses no matcher. The shared greedy registration pass
  enables joint `pose_cont` after each discrete winner.
- `post_matcher`: enable joint `pose_cont` after every main `prob_neigh`
  matcher winner.
- `standalone_final`: run the ordinary `prob_neigh` main loop, then one
  full-particle standalone `refine=pose_cont` iteration. Before that terminal
  pass, write a recoverable project checkpoint and an independent copy of its
  registered canonical sigma state. Iterative references are retained.

All modes force `inpl_cont=no` and `pose_cont_route=joint`. The public UI does
not expose `refine`, `pose_cont`, `inpl_cont`, or `pose_cont_route`; those are
child-stage controls owned by the commander.

## Shared lifecycle

The workflow reuses the `refine3D_auto` startup, optional external-reference
pose initialization, masked greedy registration, reconstruction handoffs, and
native-sampling final reconstruction. Child invocations use `mkdir=no` and
monotonically increasing iteration numbers. The final reconstruction consumes
the last child stage's project and reconstruction state.

The standalone Cartesian route retains its existing Cartesian sigma
contribution and does not require live PFTC state. Shared-memory and
distributed execution use the same commander sequence and project handoffs.

Related policies:

- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D_policy.md](refine3D_policy.md)
- [Standalone Cartesian implementation note](../../implementation_notes/continuous_3D_pose_cont_refine3D.md)
