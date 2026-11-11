# Local 3D Refinement

## Problem

A discrete 3D refinement ([refine3D](refine3d.md), run by [solve3D](solve3d.md)
or `refine3D_auto`) leaves every particle with a pose: a projection direction
from a quasi-uniform grid, an in-plane angle and a shift. Its precision is
bounded by the grid. `nspace = 20000` spaces directions about 1.4 degrees
apart, and the polar references fix the in-plane sampling. At high resolution
the pose error has to be well below that spacing.

`refine3D_local` continues such a refinement without a grid. For one state,
from the stored poses and the existing half maps, it refines each particle's
five pose parameters continuously and locally:

- three rotation components;
- two shift components.

It does not search globally, and it does not assign states. The output is a
refined set of poses and a map reconstructed from them at native sampling.

## Model

**Pose.** A pose is a rotation `R` and a shift `s`. A rotation is updated on
the right, `R <- R exp([w]x)`, where `w` is a 3-vector and `[w]x` its skew
matrix. The parametrization has no singularity near the current pose: unlike
Euler angles, every direction of rotation costs the same.

**Model of one particle.** By the projection-slice theorem, a particle's
Fourier pixel `(h, k)` samples the 3D transform of its half's map on the
rotated plane, `q = padf [h, k, 0] R exp([w]x)`. The model of the pixel is

```text
m(h,k) = |CTF(h,k)| F_half(q) exp(2 pi i (h s1 + k s2)/box)
```

with these conventions:

- The particle is phase-flipped, so the CTF enters as its modulus.
- The particle is pre-shifted by its stored shift, so `s` is an increment on it.
- The reference is the half map of the particle's own half (gold standard),
  padded twofold.

**Objective.** The Euclidean objective of the series, normalized by the
particle's own whitened power:

```text
L(w, s) = sum_(h,k) |y(h,k) - m(h,k)|^2 / sigma2(k)  /  sum_(h,k) |y(h,k)|^2 / sigma2(k),
```

The sums run over the Cartesian pixels of the current band, and `sigma2(k)` is
the noise power of the particle's half at shell `k`. Cartesian pixels already
count each shell in proportion to its circumference, so the `k` weight of the
polar form is implicit. The score `exp(-L)` is stored per particle. It is not
comparable with the polar score.

The steps below describe the default `cont_prior=none` case. Optional
`cont_prior=shift|pose` adds a seed-centred Gaussian shift penalty and, for
`pose`, a Matrix-Fisher rotation penalty. Acceptance and stopping use the
combined objective; the reported score and noise residual remain data-only.
The fixed shift width `cont_prior_shift_sigma` is a per-axis standard deviation
in Angstroms, matching the polar command-line convention; zero selects the
automatic estimator. Internally, the prior uses native-pixel variance.
The [delivered contract](../implementation_notes/completed/refine3d_local_cartesian_priors_done.md#delivered-implementation)
and [open concerns and plan](../implementation_notes/planned/refine3d_local_cartesian_priors_todo.md)
describe the current widths and proposed posterior-derived search restraints.

**Differentiable sampling.** The model must be differentiable in the pose,
including across grid cells. `F` is therefore sampled off grid with a
separable quadratic B-spline, three taps per axis. Its value and first
derivative are continuous at every cell switch. The truncated Kaiser-Bessel
(KB) stencil of the polar searches is not: its value jumps by about 1 % of the
outer coefficients' difference where a tap leaves the window.

Two choices keep the Cartesian and polar representations consistent, so both
share one noise model:

- The reference coefficients are weighted by the ratio of the KB stencil's
  response to the B-spline's. On the grid, the B-spline gather then reproduces
  the KB central sections of the polar searches.
- The particle carries the KB stencil's real-space taper, as the polar
  particles do.

## Algorithm

### One particle: a bounded Levenberg-Marquardt transaction

Each sampled particle runs one transaction from its stored pose `(R_0, s_0)`.

1. **Evaluate.** One pass over the particle's band returns `L`, its gradient
   `g` in the five parameters, and the Gauss-Newton matrix `H = J^T J` (scaled
   as for `L`).
   - The rotation columns of the Jacobian `J` come from the B-spline gather's
     three spatial derivatives, by the chain rule through `q`.
   - The shift columns come from the phase ramp.
   - For a real image every summand is Hermitian between `(h, k)` and
     `(-h, -k)`, so the pass runs over half the plane and doubles the sum.
2. **Scale the curvature.** `H` is multiplied by the regression slope of the
   data on the model, weighted by `|f|^2` as `H` is, and clamped to
   `[0.05, 20]`. A model brighter than the data would otherwise overstate the
   curvature and the step would overshoot.
3. **Propose.**
   - The rotation coordinates are scaled by 0.1 rad.
   - The step solves `(H_s + mu D) d = -g_s` by Cholesky, where `D` is the
     diagonal of the scaled matrix floored relative to its largest element,
     and `mu` starts at `1e-3`.
   - Each proposal is capped at 0.1 rad of rotation and one native pixel of
     shift.
   - Over the whole transaction, the pose may move at most `athres_cont`
     (default 10 degrees) of rotation and `trs` of shift from `(R_0, s_0)`.
     A proposal beyond either bound is rejected without evaluation.
4. **Accept or damp.** With the predicted reduction
   `p = -g^T d - (1/2) d^T H d`, a trial pose is evaluated and accepted when
   the actual reduction is positive and at least `0.25 p`.
   - After an acceptance, `mu` halves when the ratio exceeds 0.75.
   - After a rejection, `mu` grows by a factor that starts at 4 and doubles
     with each consecutive rejection.
5. **Stop.** The transaction ends:
   - when an accepted step reduces `L` by less than `1e-6` of `L`;
   - when a proposal predicts less than that (it is not evaluated);
   - after 8 consecutive rejections;
   - after 40 iterations;
   - when the gradient vanishes, or the system is singular or ill-conditioned.
6. **Commit.** Only a finite pose with a lower objective is written. Any other
   outcome returns the stored pose unchanged, with a typed status (no
   improvement, unreliable system, bound-rejected, non-finite).

The per-shell residual of the evaluation that produced the committed pose is
the particle's contribution to the next noise spectrum, at no extra pass. On
the beta-galactosidase benchmark a solve takes 4.7 passes on average.

### One iteration

The iteration is the [refine3D](refine3d.md) round, with the transaction above
as its search:

1. **Sample.** Choose the particles to update. Particles that have been
   updated least are drawn first.
2. **Prepare references.** Use the even and odd maps of the last assembly: NU
   filtered, ML regularized, and masked as in refine3D. The matching band runs
   to the nonuniform-filtering handoff: the finest member of the NU bank
   retained under its FSC-derived cap. `lpstop` and the crop's Nyquist cap the
   band further.
3. **Refine and commit.** Run one transaction per sampled particle against the
   map of its own half.
4. **Update the noise model** from the committed residuals.
5. **Reconstruct.** Insert the sampled particles at their new poses, blend
   them with the trailing reconstruction of the earlier iterations, and
   assemble.
6. **Assemble.** Compute the FSC, apply ML regularization, NU filtering and
   automasking, and set the next band.

After the last iteration, all particles are reconstructed at native sampling,
and the map is postprocessed.

### Starting from a finished refinement

The command builds nothing from the particles before the first iteration.

It requires:
- the merged map and its even/odd pair;
- stored poses for every active particle, with half labels, and particles in
  both halves.

It takes over, when compatible:
- the noise spectra;
- the NU products;
- the trailing reconstruction.

A missing FSC is computed from the halves. Without a compatible trailing
reconstruction, the first iteration updates every particle and seeds one; later
iterations update the sampled fraction.

### Schedule and fixed choices

| Setting | Value |
| --- | --- |
| Objective | Euclidean, with ML regularization |
| Halves | kept separate, each particle against its own half |
| Reference centring | none |
| Sampling | native by default (`autoscale=no`) |
| Particles per iteration | `nsample` (default 25000), as the fraction `nsample / N` of the `N` active particles; every particle when that fraction exceeds 0.99 |
| Iteration budget | about four updates per particle, `4 x ceil(N / nsample)` iterations, between 3 and 50 |
| Shift bound `trs` | `clamp(12 A / smpd, 2, 8)` pixels unless given |

The command refuses a projection grid, a fixed `lp`, the LP-set NU mode, more
than one state, and any setting that changes the Euclidean scoring or merges
the halves.

### A polar local search, for comparison

With `refine=prob_neigh` the command keeps everything above but the search. In
place of the transaction, the polar machinery of [refine3D](refine3d.md) runs
in its `geom` mode: it scores the projection directions within `athres_cont`
(the same 10 degrees) of each stored pose, every in-plane rotation of each,
then refines the shift under the shift prior and refines the in-plane angle
continuously, and commits a draw from the tempered posterior over the scored
candidates. Directions stay on the grid (`nspace = 20000`, about 1.4 degrees).
The probabilistic modes have no stable-particle test, so the run ends on its
iteration budget or the resolution exit; equal `maxits` makes the two searches
comparable. `lpstop` caps the band without merging them. For a multi-state project,
`state=X` refines that state on a single-state copy.

## Convergence

The run converges when consecutive passing iterations have together sampled
90 % of the active particles. A failing iteration resets the streak.

An iteration passes when, in every populated state, 90 % of the sampled
particles are stable. A particle is stable when:

- its transaction ended accepted, or with a finite objective and no
  improvement;
- it moved by at most 0.5 degrees of rotation and 1 A of shift.

A transaction that rejected every proposal keeps its stored pose. A zero pose update
alone therefore never counts as convergence. The three thresholds are
`cont_conv_rot`, `cont_conv_shift` and `cont_conv_frac`.

## Rationale

- **Continuous rather than finer.** A finer grid costs work in proportion to
  the number of directions and stays discrete. A Gauss-Newton step from a pose
  near the optimum reaches it in a few evaluations of one particle, with no
  grid at all.
- **Local by construction.** Started from a converged discrete refinement, the
  poses are near their optimum, and a global search mostly re-finds them. The
  10-degree capture range absorbs small misassignments. It does not repair
  large ones, which remain the discrete search's job.
- **Five parameters jointly.** The rotation and shift are coupled through the
  projection. The joint stage alone needs 37 % fewer iterations, and 38 % less
  time, than a shift-only stage followed by the joint one.
- **Commit only improvements.** Each particle's objective cannot rise within
  an iteration. The typed outcome separates a converged particle from an
  unobservable one or one held at a bound, and convergence counts only the
  first kind.
- **One noise model.** The B-spline weighting and the KB taper make the
  Cartesian and polar objectives two samplings of one model. A local
  refinement can therefore start from the noise spectra of the discrete one,
  and hand its own back.
- **Measured cost.** On the beta-galactosidase benchmark, the orientation
  search fell from 19.3 s to 1.9 s at 16 threads. Two kinds of change made
  the gain:
  - **Fewer passes per solve, 13.75 down to 4.7:** duplicate passes removed,
    the curvature scale, the predicted-reduction stop, and the residual taken
    from the accepted evaluation.
  - **Cheaper passes, 215 down to 62 ns per pixel evaluation:** the continuous
    gather, which also lowered the rotation error, the half plane, and a
    separable gather.

## Implementation

- Command and defaults: `src/main/commanders/simple/simple_commanders_refine3D.f90`
  (`exec_refine3D_local`, `exec_refine3D_workflow`); startup preflight and
  staging: `src/main/simple_refine3D_workflow_controller.f90`.
- Iteration strategy: `src/main/strategies/parallelization/simple_refine3D_strategy.f90`.
- Per-particle search: `src/main/strategies/search/simple_strategy3D_cont.f90`.
- Objective, gradient and Gauss-Newton terms: `src/main/cftc/simple_cartft_calc.f90`;
  the Levenberg-Marquardt transaction: `src/main/cftc/simple_cartft_pose_opt.f90`.
- Convergence: `check_cont_conv` in `src/main/simple_convergence.f90`.
- Policies: `doc/policies/3D/refine3D_local_policy.md` and the continuous
  pose section of `doc/policies/3D/refine3D_policy.md`; the measured gains in
  `doc/implementation_notes/completed/cont_lm_acceleration.md`.
