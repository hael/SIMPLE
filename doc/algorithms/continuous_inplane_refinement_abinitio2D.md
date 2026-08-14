# Continuous In-Plane Refinement in `abinitio2D`

## Public interface

```text
inpl_cont=yes|no
```

`no` is the default and preserves the legacy discrete shift/angle search.
`yes` follows the polish-only principle: continuous refinement changes pose
precision, never search behavior. All selection machinery -- candidate
scoring, probabilistic table construction, shift-seed estimation -- runs the
legacy discrete route unchanged, and the joint continuous `(sx,sy,theta)`
optimizer runs exactly once per particle per iteration to polish the
committed assignment. The search trajectory is therefore identical to the
validated legacy trajectory by construction, on every sample.

Continuous refinement is active for the raw Euclidean and cc objectives (the
cc route minimizes `-cc` with a quotient-rule angular derivative and reports
the clamped correlation as its score). It is unsupported for hybrid/denoised
objectives; those opt-in requests fail rather than falling back to the
callback. `abinitio2D_sgd` rejects `inpl_cont=yes`,
while time-series shift-only search invokes neither an angle callback nor the
joint angle optimizer. The `abinitio2D` controller propagates the option to
every `cluster2D` child, including terminal probabilistic and dense passes.

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

With eligible `inpl_cont=yes`, the search runs exactly as `inpl_cont=no`
until a class and integer in-plane cell have been committed. One joint solve
then polishes the pose:

1. the selected cell is authoritative: the solve seeds at
   `(sx_seed,sy_seed,theta_selected)` with `theta` limited to plus or minus
   two grid units, and never rescans other cells;
2. one three-variable L-BFGS-B solve refines the pose against the raw
   Euclidean objective, with all three analytic derivatives evaluated at the
   same continuous pose;
3. a finite, materially improving result commits the fractional pose;
   otherwise the selected discrete pose stands, re-scored at its own shift.

Continuous initialization never performs the legacy 5-by-5
shift/all-rotation coarse scan and never seeds from fractional `e3`
metadata.

Probability tables are pure legacy under both `inpl_cont` values: candidate
scoring, shift-seed estimation, and per-candidate shift refinement use the
legacy discrete optimizer, and the table stores only the integer `inpl`,
score, and shift in the rounded-index frame. Polishing candidates would
sharpen the score contrast the stochastic sampler sees and collapse
exploration on difficult samples.

After the final probabilistic class/in-plane decision, the matcher seeds the
joint solve at the assigned integer cell and its recovered native shift.
Only this final assignment persists the accepted continuous angle in `e3`;
the nearest integer remains in `inpl` for compatibility, and the final shift
and score are written from the same joint pose. A valid non-improving solve
retains the assigned discrete pose with its re-scored objective value; a
numerically invalid solve leaves the incoming assignment untouched.

Statistics parity: `dist_inpl` and every other statistic that steers search
control or convergence is computed from the discrete cells on both sides of
the comparison; the fractional `e3` is carried for reconstruction and
persistence only. Sub-grid precision would otherwise read as premature
convergence and truncate the annealing schedule difficult samples require.

The final `store_solution` of the joint pose bypasses the score-improvement
guard: the joint result restarts from the native shift, so its score is not
comparable to the stored peak-search score (obtained at a shift that is
discarded), and rejecting the store would couple the committed joint shift
with a stale integer angle in `assign_ori`. When `assign_ori`'s class argmax
selects a class other than the jointly refined one, that class's grid angle
is written (the coupled-pose guard), preserving legacy selection semantics.
Evaluation-only passes (`refine=eval`) restore the incoming `e3` verbatim; a
pass that performs no search must not degrade a fractional angle to its grid
cell.

## Joint evaluator

The evaluator uses thread-local coefficient workspaces. For each pose it
fills a `(nrots, 3·nk)` buffer with three column sections on the polar
samples — `S·REF`, `argtransf_x·S·REF`, and `argtransf_y·S·REF`, the latter
two extended to the second half-circle as `-conjg(...)` because the phase
arguments flip sign across the Friedel mate — and forward-transforms all
three sections in a single batched FFT execution (`plan_fwd3_many`). The
resulting coefficient series for the objective and both shift gradients
evaluate directly at `theta`, and the angle gradient follows from
differentiating the series. It performs no inverse FFT, constructs no
discrete loss vector, and makes no per-evaluation allocation. At exactly
zero shift only the shift-phase multiplication is skipped; there is no
dead-zone threshold, so the evaluator is smooth in the shift variables. The `argtransf` products must be formed before the angular
transform; multiplying after it is not equivalent (it would require a
convolution in angular frequency) and silently corrupts the shift gradient.

## Shift-search API

`pftc_shsrch_grad` exposes separate constructors for its valid configurations:

- `new_legacy`: legacy L-BFGS-B shift search with the discrete angle callback;
- `new_fixed`: fixed-angle L-BFGS-B shift search;
- `new_direct`: fixed-angle bounded direct-gradient shift updates; and
- `new_joint`: joint `(sx,sy,theta)` L-BFGS-B refinement.

The callback is attached only by `new_legacy`. This makes route audits
mechanical and avoids optional mode-flag combinations. The obsolete
`shbarrier` option was also removed because it did not affect the
implementation.

`select_best_discrete_angle` is shared by the legacy callback and joint
initialization. Joint code calls it once at the supplied shift; it does not
attach the callback to L-BFGS-B.

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
