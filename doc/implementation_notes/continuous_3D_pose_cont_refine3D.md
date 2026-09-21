# Standalone Cartesian pose refinement in `refine3D`

**Contract status:** IMPLEMENTED; matcher consolidation and the compact matched-
projector regression passed compilation and runtime validation.

This living note keeps two deliberately separate routes:

- `pose_cont=yes` preserves the post-matcher transaction introduced by commit
  `65624b924`; it polishes the authoritative PFTC or `inpl_cont` winner.
- `refine=pose_cont` selects the bona fide standalone local-search strategy. It
  starts from initialized `ptcl3D` poses and owns particle-pose search without
  executing PFTC or `inpl_cont`.

The selectors are mutually exclusive. The standalone implementation is
feature-complete. The three final source-review findings are corrected, and
the affected compact and Server validation has passed.

## Scientific contract

1. `strategy3D_pose_cont` is a bona fide strategy class in
   `simple_strategy3D_pose_cont.f90` and implements the standard `strategy3D`
   lifecycle: `new`, `srch`, `oris_assign`, and `kill`.
2. The strategy is a local continuous polisher, not an ab-initio orientation
   finder. Every active particle must have a positive state, a valid half-set,
   finite Euler angles and shifts, and a positive stored projection index.
   Identity Euler angles are valid.
3. `refine3D` owns the common plumbing: reference masking, filtering and
   centering; raw particle reads; noise normalization and masking; Fourier
   cropping; CTF metadata; shell-noise weights; half-set identity; and
   reconstruction.
4. One immutable `pose_cont_reference_workspace` is shared by a matcher pass.
   Each strategy instance binds that workspace, one thread-local Cartesian
   image, and its batch-local particle index.
5. Search ownership is Cartesian only. The standalone strategy does not
   initialize or call `strategy3D_srch`, PFTC projection banks, PFTC scoring,
   discrete in-plane search, or `inpl_cont`.
6. State and half-set are fixed. The selected LM route refines the stored
   rotation and shift within its local bounds.
7. An accepted improvement is converted to SIMPLE Euler angles and
   native-pixel shifts. Invalid, bounded-out, or non-improving transactions
   preserve the seed. The Cartesian objective never replaces the PFTC `corr`.
8. Reconstruction remains unchanged and consumes the resulting project pose.
9. The production Cartesian LM objective is normalized correlation
   (`POSE_CONT_OBJECTIVE_CART_NCC`). Cartesian Euclidean remains an internal
   validation option. SIMPLE's surrounding `objfun=euclid|cc` contract is
   independent and does not select the Cartesian objective.

## Implemented ownership map

- `simple_ui_refine3D.f90` exposes the standalone `refine=pose_cont` choice and
  documents its separation from `pose_cont=yes`.
- `simple_matcher_refvol_utils.f90` reuses ordinary reference preparation but
  skips PFTC bank generation for the Cartesian-only route.
- `simple_matcher_ptcl_batch.f90` reuses ordinary particle I/O while skipping
  polarization and unused padded per-thread images.
- `simple_strategy3D_matcher.f90` validates seeds, routes the standalone class,
  binds matcher-owned context, and preserves ordinary reconstruction.
- `simple_strategy3D_pose_cont.f90` owns the standalone strategy lifecycle and
  project-pose transaction.
- `simple_pose_cont_refine3D_adapter.f90` owns reference workspaces, particle
  preparation, route selection, rollback, and Cartesian sigma contributions.
- `simple_cartesian_pose_refiner.f90` remains the numerical owner of the
  two-parameter shift and five-parameter Cartesian LM solves.
- `simple_euclid_sigma2.f90` accepts per-particle Cartesian residual-shell
  contributions without requiring a PFTC object.

## Class lifecycle

### `new`

- Retain particle identity and validate the standalone strategy contract.
- Require `oritype=ptcl3D`, `objfun=euclid`, and `inpl_cont=no`.
- Do not initialize inherited PFTC search state.

### `bind_context`

- Bind the matcher-owned reference workspace and builder/parameter context.
- Bind a thread-local cropped image and record the batch-local particle index.
- Copy the matcher-owned LM route and bounds so execution and reporting use the
  same canonical adapter policy.

### `srch`

1. Capture the stored pose as both seed and rollback point.
2. Prepare the already-loaded particle on the Cartesian Fourier grid.
3. Select shell weights and the matching state/half reference.
4. Run `shift_then_joint` or `joint` LM from the stored pose.
5. Stage an accepted pose; otherwise retain the seed.
6. Evaluate the terminal pose for sigma accounting and assign the result.

### `oris_assign`

- Commit Euler angles and native-pixel shifts only after accepted improvement.
- Preserve `corr`, state, half-set, and discrete seed fields.
- Refresh convergence telemetry from the actual seed-to-terminal motion; do
  not retain values from the project that supplied the initial pose.

### `kill`

- Release local orientation state and nullify borrowed pointers.
- Do not destroy matcher-owned images, builder state, or reference workspaces.

## Completed pre-commit corrections

1. **Compile blocker:** removed the accidental leading `n` before the
   `result%bound_hits` assignment in `add_stage_accounting`.
2. **Preserve the delivered post-PFTC contract:** `pose_cont=yes` supports
   ordinary and probabilistic pose-search modes only after a valid assignment.
   It rejects `sigma`, `eval`, and fixed-orientation `prob_state` modes rather
   than polishing a stale or intentionally fixed seed.
3. **Make standalone convergence metadata authoritative:**
   `strategy3D_pose_cont` now refreshes `dist`, `dist_inpl`, `shincarg`,
   `mi_proj`, `mi_state`, and the local-search `frac` policy for every terminal
   transaction, including rollback. The matcher writes `frac_greedy=0` for
   the standalone route. The earlier Server output inherited stale values from
   the input project, so its convergence declaration and search-statistics
   lines remain non-acceptance evidence until the smoke test is repeated.

The compact class test currently proves the identity-seed validation contract;
the complete class lifecycle is covered only by the production smoke run. A
focused lifecycle test remains desirable but is not required to diagnose the
three corrections above.

## Observed Server evidence

Evidence directory:
`pose_cont_strategy_20260916_125215`.

- The compact numerical and adapter suites passed.
- Standalone `refine=pose_cont` processed all 3,081 particles, reconstructed
  normally, and ended with `SIMPLE_REFINE3D NORMAL STOP`. Matcher time was
  595.59 seconds and peak RSS was 2.532 GiB.
- The retained `refine=shc pose_cont=yes` route and the default
  `refine=shc pose_cont=no` route also ended normally.

These three runs are functional smoke evidence, not a matched scientific or
performance comparison. The standalone run reported starting-reference
foreground sigma near `1.14e-2`, while both SHC runs reported about `1.32e-2`.
The standalone and post-PFTC runs reused canonical sigma state, whereas the
default SHC run rebuilt it. Therefore their FSC and runtime differences cannot
be attributed solely to the selected search route.

The post-PFTC handoff investigation found that PFTC preparation mutates its
particle image. The adapter now preserves the raw particle before PFTC and owns
the copy through Cartesian preparation. Matched production runs reproduced the
cropped and padded preparation boundaries exactly, and the adapter contract
tests cover Euler/matrix conversion, native/cropped shifts, accepted assignment,
rollback, and metadata preservation.

Final corrected evidence directory:
`pose_cont_strategy_20260916_150959`.

- Both compact suites passed: numerical owner, solver, complete numerical
  suite, and refine3D adapter.
- Standalone `refine=pose_cont` processed all 3,081 particles, rebuilt the
  canonical sigma state, reconstructed normally, and ended with
  `SIMPLE_REFINE3D NORMAL STOP`.
- Corrected telemetry reports zero greedy searches, 100% local-search
  coverage, 2.855 degrees mean seed-to-terminal orientation motion, 2.374
  degrees mean in-plane motion, and 0.396 pixels mean shift motion. The run is
  correctly not declared converged after this one iteration.
- Matching took 747.38 seconds, peak matcher RSS was 2.682 GiB, and the final
  FSC resolutions were 3.617 A at 0.143 and 4.322 A at 0.5.

## Matcher consolidation and compact boundary regression

The final cleanup keeps `simple_strategy3D_matcher` orchestration-only:

1. `simple_pose_cont_refine3D_adapter` constructs the canonical route policy
   and native-to-cropped shift limits for both pose-cont entry routes.
2. The adapter validates global sigma-shell and particle bounds, constructs the
   zero-based particle sigma vector, and prepares the Cartesian particle data.
   A compact `pose_cont_particle_spec` groups state, half-set, particle index,
   shell range, and CTF metadata; the large particle and sigma arrays remain
   borrowed inputs. The matcher and standalone strategy no longer duplicate
   this numerical plumbing.
3. The compact numerical suite compares the Cartesian gather with SIMPLE's
   established oversampled projector kernel at identical rotated coordinates
   from the same physical reference. Correlation, fitted gain, and relative-L2
   gates protect the matched reference boundary without production diagnostics.

Compilation and both compact suites passed on 2026-09-21. The matched-projector
regression reported correlation `1.0000E+00`, fitted gain `1.0000E+00`, and
relative L2 `1.4300E-07`; the numerical, solver, mother, and adapter suites all
reported `PASS`.

## Remaining validation and delivery order

1. Run current-`master` one-iteration smokes for default-off, post-matcher
   `pose_cont=yes`, standalone `refine=pose_cont`, and probabilistic polishing.
2. Design the separate developer `refine3D_auto_pose_cont` policy, including
   standalone particle sampling, without changing the established
   `refine3D_auto` UI.
3. Validate the final objective and acceptance policy on a harder real dataset.
