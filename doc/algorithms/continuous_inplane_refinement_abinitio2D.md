# Continuous In-Plane Refinement in `abinitio2D`

## Vocabulary

This development uses three algorithm names:

- **legacy discrete**: the historical shift optimizer and discrete angle update;
- **continuous callback**: alternating shift optimization and continuous angle updates;
- **continuous joint**: simultaneous optimization of `(sx,sy,theta)`.

“Callback” describes the alternating continuous algorithm because its angle
update is executed through the shift optimizer's callback. The historical
algorithm is called *legacy discrete*, even though its discrete angle update is
also implemented using callback plumbing.

Candidate selection and refinement are separate concepts. Candidate selection
remains discrete for every route. The two active routes refine the selected
candidate continuously.

## Public interface

```text
inpl_cont=no|callback|joint
```

| Value | Algorithm | Final angle |
|---|---|---|
| `no` | Legacy discrete shift/angle search | Discrete grid angle |
| `callback` | Continuous callback refinement | Accepted continuous angle |
| `joint` | Continuous joint refinement | Accepted continuous angle |

The routes are mutually exclusive. `callback` never invokes the joint
optimizer, and `joint` never invokes a continuous angle callback.

`no` is the default. It constructs neither continuous implementation and
preserves the historical algorithm.

## Eligibility boundary

Continuous refinement is available only when all of the following are true:

- `objfun=euclid`;
- `.not. l_objfun_den`;
- streaming SGD is inactive; and
- time-series shift-only search is inactive.

When an active mode is requested outside this boundary, the effective mode is
`no`. Low-level constructors and evaluation gateways also reject continuous
use without the raw Euclidean objective. A continuous derivative for the
hybrid/denoised objective has not been implemented.

The `abinitio2D` controller copies `inpl_cont` to every `cluster2D` child
command, including the final probabilistic and dense all-particle passes.

## Shared discrete candidate selection

All routes first use the existing polar machinery to select a class and an
integer in-plane rotation. Probabilistic stages still construct their
probability table and draw an integer class/rotation assignment exactly as
before.

For a controlled comparison, both continuous routes then start from:

- the same selected class;
- the same selected integer rotation;
- the same native-frame shift seed; and
- the same ordinary shift limits.

When a previous-state shift seed is used, its search is fixed at the previous
discrete angle for both continuous routes. It cannot perform an earlier angle
update that would make the starting states different.

## `inpl_cont=no`: legacy discrete

For a selected particle/reference pair, the legacy route is:

1. Evaluate the polar objective over integer rotation indices.
2. Select the best discrete angle.
3. Run two-variable L-BFGS-B over `(sx,sy)` at that integer angle.
4. Re-evaluate all integer angles and update the selected index.
5. Alternate shift and angle updates for at most five cycles, stopping when the
   integer angle is unchanged.
6. Store the integer index in `inpl` and reconstruct `e3` from that index.

No continuous state, continuous callback, or joint optimizer participates in
this route.

## `inpl_cont=callback`: continuous callback

The callback route retains the alternating structure but makes the angle
coordinate continuous:

1. Start at the selected integer rotation and native-frame shift seed.
2. Build the Euclidean angular Fourier coefficients at the current shift.
3. Starting from the selected angle, perform up to three safeguarded Newton
   updates in the continuous angular coordinate.
4. Run two-variable L-BFGS-B over `(sx,sy)`.
5. Invoke the continuous angle update through the optimizer callback and
   alternate until the rounded angle index is unchanged or the cycle limit is
   reached.
6. Accept the result when it improves the selected discrete pose. A pure angle
   improvement is sufficient; a shift change is not required.
7. Rotate the shift with the accepted continuous angle and store the continuous
   value in `e3`.

The two shift derivatives are evaluated at the nearest integer angle. Thus the
callback route is an alternating continuous-angle method, not a coupled
three-variable derivative calculation.

## `inpl_cont=joint`: continuous joint

The joint route removes the continuous angle callback:

1. Start at `(sx_seed,sy_seed,theta_grid)` for the same selected candidate.
2. Set the angular limits to plus or minus two grid units around that candidate.
3. Run one three-variable L-BFGS-B solve over `(sx,sy,theta)`.
4. Evaluate the raw Euclidean loss and all three analytic partial derivatives
   at the same continuous pose.
5. Accept only a finite result that improves the selected discrete pose by more
   than a round-off-scaled tolerance.
6. Keep the nearest integer index in `inpl` for compatibility and store the
   accepted continuous angle in `e3`.

This is the coupled formulation: shifts and angle are components of one
optimization variable and one gradient.

## Probabilistic selection and final class averaging

Probabilistic selection remains discrete and unchanged. After it selects one
candidate, the selected continuous route is applied only to that candidate:

- `callback` runs continuous callback refinement;
- `joint` runs continuous joint refinement; and
- `no` performs no additional refinement.

The selected reference-frame shift is converted back to the native polar frame
before refinement. On acceptance, `inpl` retains the nearest integer index and
`e3` retains the continuous angle.

Orientation assignment completes before online class-average restoration reads
the particle metadata. The class averager reads real-valued `e3` directly and
passes it to `rotmat2d`; it does not reconstruct the angle from `inpl`. The same
metadata is consumed by the terminal class-average generation. Consequently an
accepted continuous angle survives through the final `abinitio2D` products.

## Continuous accuracy and polar efficiency

### Representation efficiency

Both continuous methods retain the important polar representation:

- discrete candidate discovery remains a fast polar scan;
- no Cartesian particle image is rotated during refinement;
- continuous angles are evaluated from the band-limited angular Fourier
  series;
- angular periodicity follows directly from that series;
- analytic derivatives are used; and
- continuous refinement is local to one already-selected candidate.

Continuous accuracy therefore does not require a finer global rotation grid or
repeated Cartesian interpolation.

### Callback evaluator cost

One continuous callback update builds the angular coefficient series at the
current shift and evaluates that series at a few Newton trial angles. The
coefficient construction dominates; subsequent angle trials are inexpensive.
The shift optimizer remains two-dimensional but may request several callback
updates during one alternating search.

### Current joint evaluator cost

The current joint objective/gradient evaluator is correct but not yet the most
efficient polar implementation. Each evaluation currently:

1. allocates one `nrots` loss vector and three coefficient arrays;
2. constructs coefficients separately for the objective, x derivative, and y
   derivative;
3. at nonzero shift, repeats shifted-reference preparation and the forward
   polar-ring FFT for all three constructions;
4. performs three inverse FFTs to form discrete loss vectors that the joint
   optimizer does not consume; and
5. discards the three loss vectors after evaluating their coefficient series at
   one continuous angle.

This can amount to three shifted-reference forward FFT batches, three inverse
FFTs, three coefficient passes, and four allocations per joint evaluation.

### Intended coefficient-only joint evaluator

The same objective can be evaluated more efficiently by:

- reusing thread-local coefficient workspaces;
- preparing and transforming the shifted reference once;
- accumulating objective, x-derivative, and y-derivative coefficients in one
  polar-ring traversal;
- skipping all inverse FFTs and discrete loss vectors; and
- evaluating the three coefficient series directly at `theta`.

At zero shift, the memoized reference transform can be reused. At nonzero
shift, the target is one shifted-reference transform, no inverse FFT, and no
per-evaluation allocation.

Current timing therefore compares the implemented algorithms, but it is an
upper bound on the production cost of the continuous joint formulation.

## Validation

Focused tests cover:

- default low-level construction remains discrete;
- route identity for `no`, `callback`, and `joint`;
- `callback` constructs the continuous callback but not the joint optimizer;
- `joint` constructs the joint optimizer but not the continuous callback;
- both continuous routes use fixed-angle previous-state shift seeding;
- probabilistic mode retains the selected continuous implementation;
- hybrid/denoised Euclidean requests fall back to effective `no`;
- synthetic recovery reports continuous callback and continuous joint results
  independently;
- the joint analytic gradient is compared with central differences; and
- final metadata must remain grid-aligned for `no` and contain continuous
  angles for both active modes.

The data-dependent executables require external MRC/SIMPLE fixtures. A focused
build establishes interface and compilation consistency, not numerical legacy
identity or production throughput.

## Remaining issues

1. Persisted boundary angles are not canonicalized to one explicit degree
   interval.
2. The callback Newton step limit is per iteration, allowing cumulative motion
   beyond half a grid interval from its initial angle.
3. The joint evaluator needs the coefficient-only optimization described above
   before final performance conclusions.
4. Paired fixed-seed `abinitio2D` runs are still needed to verify numerical
   identity for `no` and compare accuracy, convergence, final class averages,
   and throughput for `callback` and `joint`.

## Primary code locations

- Public parameter: `src/main/params/simple_parameters*.f90` and
  `src/main/ui/simple/simple_ui_cluster2D.f90`.
- Stage propagation: `src/main/simple_abinitio2D_controller.f90`.
- Route selection and continuous pose persistence:
  `src/main/strategies/search/simple_strategy2D_srch.f90`.
- Probabilistic post-selection refinement:
  `src/main/strategies/search/simple_strategy2D_prob.f90`.
- Callback and joint optimizers: `src/main/pftc/simple_pftc_shsrch_grad.f90`.
- Continuous residual and joint gradient:
  `src/main/pftc/simple_polarft_corr.f90`.
- Focused tests: `production/tests/simple_test_euclid_route_identity.f90`,
  `production/tests/simple_test_euclid_stage1_validation.f90`, and
  `production/tests/simple_test_euclid_2d_metadata.f90`.
