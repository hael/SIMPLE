# Continuous In-Plane Refinement in `abinitio2D`

## Public interface

```text
inpl_cont=yes|no
```

`no` is the default and preserves the legacy discrete shift/angle search.
`yes` enables joint continuous refinement of `(sx,sy,theta)` after the
discrete polar search selects a candidate. There is no alternating continuous
callback route.

Continuous refinement is active only for the raw Euclidean objective. It is
disabled for hybrid/denoised objectives, streaming SGD, and time-series
shift-only search. The `abinitio2D` controller propagates the option to every
`cluster2D` child, including terminal probabilistic and dense passes.

## Discrete mode

With `inpl_cont=no`, the historical search:

1. evaluates the polar objective over integer rotation indices;
2. optimizes `(sx,sy)` at the selected integer angle;
3. re-evaluates the integer angle grid; and
4. alternates those updates until the index is unchanged or the cycle limit is
   reached.

The integer index is stored in `inpl`, and `e3` is reconstructed from that
index. No joint optimizer or continuous pose state is constructed.

## Continuous mode

With eligible `inpl_cont=yes`, candidate selection remains discrete. The one
selected candidate is then refined as follows:

1. start from `(sx_seed,sy_seed,theta_grid)`;
2. limit `theta` to plus or minus two grid units around the candidate;
3. run one three-variable L-BFGS-B solve;
4. evaluate the raw Euclidean objective and all three analytic derivatives at
   the same continuous pose; and
5. accept only a finite, round-off-significant improvement.

The nearest integer angle remains in `inpl` for compatibility, while the
accepted continuous angle is stored in `e3`. Probabilistic sampling therefore
remains discrete; only its selected result receives continuous refinement.

## Joint evaluator

The evaluator uses thread-local coefficient workspaces. For each pose it
prepares the shifted reference once, accumulates the objective, x-gradient,
y-gradient, and angle-gradient coefficients in one traversal, and evaluates
them directly at `theta`. It performs no inverse FFT, constructs no discrete
loss vector, and makes no per-evaluation allocation.

## Shift-search API

`pftc_shsrch_grad` exposes separate constructors for its valid configurations:

- `new`: legacy L-BFGS-B shift search with an optional discrete angle update;
- `new_direct`: fixed-angle bounded direct-gradient shift updates; and
- `new_joint`: joint `(sx,sy,theta)` L-BFGS-B refinement.

This avoids combinations of optional mode flags. The obsolete `shbarrier`
option was also removed because it did not affect the implementation.

## Primary code locations

- Parameter and UI: `src/main/params/simple_parameters*.f90` and
  `src/main/ui/simple/simple_ui_cluster2D.f90`.
- Stage propagation: `src/main/simple_abinitio2D_controller.f90`.
- Search selection and continuous metadata:
  `src/main/strategies/search/simple_strategy2D_srch.f90`.
- Joint optimizer: `src/main/pftc/simple_pftc_shsrch_grad.f90`.
- Coefficient-only objective/gradient:
  `src/main/pftc/simple_polarft_corr.f90`.
- Focused validation is orchestrated by
  `production/tests/simple_test_continuous_inplane_rotation2D.f90` and implemented in
  `production/tests/simple_test_continuous_inplane_rotation2D_route_identity.f90`,
  `production/tests/simple_test_continuous_inplane_rotation2D_stage1_validation.f90`, and
  `production/tests/simple_test_continuous_inplane_rotation2D_metadata.f90`.
