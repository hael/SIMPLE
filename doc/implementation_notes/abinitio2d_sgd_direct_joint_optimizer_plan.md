# Direct joint SGD optimizer for `abinitio2D_sgd`

## Status

Long-term global direction only. This is not the current implementation plan.
The immediate incremental plan is
`abinitio2d_sgd_seeded_lbfgsb_plan.md`: first give the established L-BFGS-B
solver the hard-stream winner, remove its redundant wrapper-level rotation
rescan, and measure that change before replacing the solver.

Repository baseline reviewed on 2026-08-10: `master` at `c769408f1`
(`more inpl_cont=yes updates`). The plan below incorporates the corrected
continuous-angle gradient and the current joint-optimizer contract.

This plan must be reviewed and approved before any Fortran source, test, UI,
or command behavior is changed. The existing continuous in-plane refinement
implementation is treated as a validated, independent feature and must remain
numerically and operationally unchanged.

### Implementation readiness (2026-08-07)

There is no known theoretical or numerical blocker in the proposed direct
method itself. Its implementation is nevertheless deferred until the
hard-stream-seeded L-BFGS-B A/B experiment establishes how much benefit comes
from removing the rescan while retaining the mature solver. The
corrected continuous-angle objective and derivatives are already integrated in
2D and have also been transferred to the `abinitio3D` and `refine3D` paths.
That work reduces implementation risk because the authoritative joint
Euclidean objective, angle convention, and durable orientation metadata are
no longer speculative.

The command alone does not complete the plan, however. Merely allowing

```text
simple_exec prg=abinitio2D_sgd objfun=euclid inpl_cont=yes
```

would still enter the established L-BFGS-B continuous-refinement machinery
unless Phases 2-4 below add and route the stream-owned direct joint optimizer.
Completion means that an SGD-active iteration uses `minimize_joint_raw`, while
SGD-inactive/warmup iterations and standard `abinitio2D` retain their current
routes, and that the focused plus Oracle/Linux validation gates pass.

When this global direction is resumed, the remaining gates are implementation
and validation work, not unresolved design choices:

1. freeze the corrected continuous baseline and seed-equivalence behavior;
2. implement the bounded direct joint optimizer without L-BFGS-B;
3. route it only from SGD-active `abinitio2D_sgd` iterations;
4. expose and test the explicit development-command activation; and
5. compile locally, run focused tests, then validate the scientific workflow
   on Oracle/Linux.

Initial values for angle step radius and stationarity/backtracking tolerances
remain empirical calibration parameters. They can be tuned from focused
diagnostics and are not a blocker to implementing the already specified
algorithm. A table-free `prob_snhc` analogue remains future work and is not a
prerequisite for v1.

## Objective

Allow the development command

```text
simple_exec prg=abinitio2D_sgd objfun=euclid inpl_cont=yes
```

to combine:

1. table-free hard class/rotation selection from raw Euclidean losses; and
2. bounded block-scaled first-order refinement of shift and continuous
   in-plane angle for the selected hard candidate.

The new optimizer is the SGD replacement for L-BFGS-B whenever the
`abinitio2D_sgd` iteration policy marks streaming SGD active. During
SGD-inactive stages, boundary warmups, and terminal conventional passes,
`abinitio2D_sgd` continues through the established optimizer selected by that
route. Standard `abinitio2D` remains unchanged. The direct joint optimizer must
reuse the established raw Euclidean objective and analytical derivatives rather
than wrap or call L-BFGS-B.

## Non-negotiable compatibility contract

The implementation is acceptable only if these route states remain distinct:

| Command and options | Required behavior |
| --- | --- |
| `abinitio2D`, default `inpl_cont=no` | Completely unchanged |
| `abinitio2D objfun=euclid inpl_cont=yes` | Existing continuous L-BFGS-B route, unchanged |
| `abinitio2D_sgd`, default `inpl_cont=no` | Existing hard stream plus shift-only `minimize_direct`, unchanged |
| `abinitio2D_sgd objfun=euclid inpl_cont=yes`, SGD inactive for the current iteration | Existing conventional continuous L-BFGS-B route |
| `abinitio2D_sgd objfun=euclid inpl_cont=yes`, SGD active for the current iteration | Hard stream plus direct joint SGD; L-BFGS-B is replaced |

Additional invariants:

- `inpl_cont` remains `no` by default everywhere.
- Standard `abinitio2D` receives no SGD-specific UI controls or commander
  logic.
- The existing `grad_shsrch_minimize_joint` continues to use L-BFGS-B for
  conventional and SGD-inactive iterations and returns its established score
  representation.
- The new direct joint routine must not allocate, construct, or invoke an
  L-BFGS-B optimizer.
- `objfun=cc`, hybrid/denoised objectives, time-series search, stages where
  SGD is inactive, and the conventional terminal probability pass remain on
  their existing routes.
- No SoftMax, top-K table, `eulprob_tab2D`, or likelihood-table machinery is
  introduced into the active streaming assignment.
- An invalid initial/current accepted evaluation preserves the complete
  incoming state. A non-finite speculative trial is rejected transactionally
  and cannot erase an earlier accepted finite improvement.

The two per-iteration joint routes intentionally have different seed ownership:

- an SGD-inactive `inpl_cont=yes` iteration owns its conventional
  initialization and performs one all-angle discrete scan at the supplied
  shift before L-BFGS-B; and
- an SGD-active `inpl_cont=yes` iteration reuses the class/rotation winner
  already selected by the table-free hard stream and must not repeat that
  global angle scan.

This difference is an orchestration choice, not a difference in the loss being
optimized. Both routes evaluate the same authoritative
`L(c,rho,sx,sy)` through `gen_raw_euclid_grad_at_angle`.

## Existing numerical components

### Discrete streaming selector

The stream currently evaluates

```fortran
call pftc%gen_raw_euclid_vals(iref, iptcl, shift, losses)
```

and selects one hard class/rotation winner:

\[
(c^*,r_0^*) = \operatorname*{argmin}_{c,r} L(c,r,\mathbf{s}_0).
\]

This selection remains unchanged.

### Existing continuous joint objective

`gen_raw_euclid_grad_at_angle` evaluates one normalized raw Euclidean loss
and its analytical gradient at a fractional polar rotation coordinate:

\[
L(c,\rho,\mathbf{s}),
\qquad
\nabla L =
\left(
\frac{\partial L}{\partial s_x},
\frac{\partial L}{\partial s_y},
\frac{\partial L}{\partial \rho}
\right).
\]

Here `rho` is a continuous rotation-grid coordinate. The objective is periodic
in rotation, so a local interval may cross rotation index `1` or `nrots`.

This PFTC routine is the authoritative mathematical API for both continuous
joint optimizers. Its formula must not be duplicated in the SGD implementation.

At the reviewed baseline, its corrected implementation has the following
numerical contract:

- shift derivatives are formed in the angular-sample domain from
  `argtransf_x*S*REF` and `argtransf_y*S*REF`;
- the derivative products use the anti-Friedel extension `-conjg` on the
  second half-circle;
- the shift fast path is gated only at exact zero. A `SHERRSQ`-style small-shift
  dead zone would make the evaluator nonsmooth and must not be reintroduced;
- the objective and all three derivatives are evaluated from coefficient
  series at fractional `rho` without an inverse FFT or per-evaluation temporary
  allocation; and
- evaluation is periodic in `rho` and agrees with
  `gen_raw_euclid_grad_for_rot_8` at integer rotation indices within numerical
  tolerance.

The repository already contains central-finite-difference coverage for all
three gradient components at both optimized and deliberately nonstationary
off-grid poses, plus integer-angle parity checks at several nonzero shifts.
Those tests are the baseline for this plan; the direct joint optimizer must
reuse this evaluator rather than introduce another derivative implementation.

### Existing L-BFGS-B joint optimizer

`grad_shsrch_minimize_joint` currently:

1. discards any prior angle state and calls `select_best_discrete_angle` once
   at the caller-supplied native-frame shift;
2. initializes `(sx,sy,rotind_frac)` from that selected grid cell;
3. sets a local angular interval of `selected_irot +/- 2` and invokes
   L-BFGS-B;
4. distinguishes numerical validity from strict improvement;
5. returns the selected discrete seed unchanged when evaluation is valid but
   no improving continuous step is found;
6. returns `irot=0` only when the joint evaluation is numerically invalid;
7. converts an improved fractional rotation into the nearest periodic integer
   rotation index; and
8. returns the established `exp(-L)` score and the requested shift frame.

That contract remains unchanged.

`select_best_discrete_angle` evaluates the existing vector objective for every
rotation at exactly one supplied shift and has focused route-identity coverage
against the legacy callback. It is part of conventional continuous-route seed
ownership; it is not the production initializer for the proposed stream-owned
joint optimizer.

### Existing shift-only direct optimizer

`grad_shsrch_minimize_direct` implements projected gradient descent with a
bounded shift box, backtracking, finite-value checks, monotonic acceptance,
and exact no-improvement preservation. It returns the stream merit `-L`.

The new joint direct optimizer should follow these safety rules but extend the
state from two to three coordinates.

## Proposed API and ownership

All numerical work remains in:

```text
src/main/pftc/simple_pftc_shsrch_grad.f90
```

Add the approved dedicated search mode and type-bound entry points:

```fortran
integer, parameter :: SHSRCH_JOINT_RAW = 4

procedure :: new_joint_raw => grad_shsrch_new_joint_raw
procedure :: minimize_joint_raw => grad_shsrch_minimize_joint_raw
procedure :: is_direct_joint
```

`grad_shsrch_minimize_joint_raw` is the stream-owned optimizer requested by
this plan. It must not call `self%opt_obj%minimize`.

`new_joint_raw` may reuse `opt_spec` only as an existing bounds/state container,
as `new_direct` currently does. It must not call `opt_factory%new`, allocate an
optimizer object, or install L-BFGS-B callbacks. `is_direct_joint` must identify
the dedicated mode so constructor/route tests can prove this separation.

Its seed contract is explicit rather than implicit:

```fortran
function grad_shsrch_minimize_joint_raw(self, irot, xy_in, &
    eta_shift, eta_inpl, max_steps, sh_rot, evaluation_valid, improved, &
    accepted_steps, objective_initial, objective_final, rotind_frac) result(cxy)
```

The exact argument order may follow local Fortran style during implementation,
but the selected integer `irot` and native-frame `xy_in` are the complete v1
pose seed. The routine initializes `rho=real(irot,dp)` and must not call
`select_best_discrete_angle` or rescan all angles. Do not add a redundant
`rotind_seed` until a real production caller needs a fractional restart seed;
if such an API is added later, its preservation and window-centering semantics
require a separate contract.

The phrase "reuse `grad_shsrch_minimize_joint`" means sharing private,
side-effect-controlled building blocks with it. The direct routine must not
call the L-BFGS-B routine as a sub-optimizer.

## Private shared refactor

Do not begin with a broad refactor of the now-corrected and tested conventional
joint route. The two optimizers share an objective but deliberately do not
share initialization policy. Extract only small, side-effect-controlled helpers
when the direct implementation would otherwise duplicate established mechanics.
Candidate responsibilities are:

```fortran
evaluate_joint_raw(...)
    call gen_raw_euclid_grad_at_angle
    update objective/gradient profiling counters
    return finite loss and gradient

project_joint_trial(...)
    independently scale the shift and angle gradient blocks
    project shift and unwrapped rho into their local bounds

finalize_joint_coordinates(...)
    map rho to the nearest periodic integer rotation index
    retain unwrapped rho only as optimizer-local state
    canonicalize the accepted periodic angle for durable e3 storage
    rotate the shift only when sh_rot is requested
```

The existing `new_joint` and `minimize_joint` need not be modified merely to
force helper reuse. If a helper is extracted from them, their observable
behavior must remain unchanged. `new_joint_raw` and `minimize_joint_raw` own a
separate block-scaled first-order loop and raw-loss return contract while
calling the same PFTC objective API.

Any modification to the conventional implementation must be protected by an
equivalence test before the new route is made callable. For fixed inputs, the
pre-refactor and post-refactor L-BFGS-B route must agree within numerical
tolerance for:

- acceptance or rejection;
- integer rotation;
- fractional rotation;
- returned shift;
- returned `exp(-L)` score.

If no conventional code is changed, this refactor-equivalence gate is not
needed; the existing continuous validation and route-identity tests remain the
regression guard.

## Direct joint optimization algorithm

The workflow remains stochastic because it operates on the iteration's sampled
particle mini-batch. The three-variable inner solve is a **bounded block-scaled
projected-gradient method with Armijo backtracking**. Its unprojected search
step solves the linearized objective over a product block trust region, but the
overall method is not a classical trust-region algorithm: it does not update
its radii from agreement between predicted and actual reduction.

The caller supplies the already-selected hard-stream state

\[
(c^*,r_0^*)=\arg\min_{c,r}L(c,r,\mathbf{s}_0).
\]

The direct joint optimizer fixes `c=c*`, starts from
`rho_0=real(r_0*)`, and does not enumerate any other class or rotation. Before
production activation, a focused seed-equivalence regression must compare, at
identical `(class,particle,shift)` inputs, the stream winner from
`gen_raw_euclid_vals` with the integer winner returned by
`select_best_discrete_angle` under raw Euclidean mode. Require the same index
when the best-versus-second-best loss margin exceeds a scaled numerical
tolerance. For a tie or near tie, require the two selected losses to be
equivalent within tolerance and document the deterministic candidate-order
rule. This validates skipping the redundant scan without making a brittle claim
that floating-point argmins are always unique.

Let

\[
\mathbf{z}=(s_x,s_y,\rho),
\qquad
\mathbf{g}=\nabla L(\mathbf{z}).
\]

The legal region is:

\[
s_x,s_y\in[-\mathrm{trs},\mathrm{trs}],
\qquad
\rho\in[r_0-w,r_0+w],
\]

where the initial proposal is `w=2` rotation bins, matching the established
continuous route. The PFTC evaluation remains periodic when this unwrapped
interval crosses a grid boundary.

### Block scaling and first-order trust region

Shift coordinates are measured in pixels; `rho` is measured in rotation bins.
They must not share an unscaled three-dimensional gradient norm.

The descent direction normalizes the two physical blocks independently. In the
interior, it is the exact minimizer of the local linear model

\[
\min_{\Delta\mathbf{s},\Delta\rho}
\mathbf{g}_s^T\Delta\mathbf{s}+g_\rho\Delta\rho
\]

subject to

\[
\lVert\Delta\mathbf{s}\rVert_2\leq\eta_s,
\qquad
|\Delta\rho|\leq\eta_\rho.
\]

The corresponding normalized direction is:

\[
\mathbf{d}_s =
\begin{cases}
\mathbf{g}_s/\lVert\mathbf{g}_s\rVert,&
\lVert\mathbf{g}_s\rVert>\epsilon,\\
0,&\text{otherwise},
\end{cases}
\]

\[
d_\rho =
\begin{cases}
\operatorname{sign}(g_\rho),&|g_\rho|>\epsilon,\\
0,&\text{otherwise}.
\end{cases}
\]

A trial with common backtracking factor `alpha` is

\[
\mathbf{s}_{trial}
=\Pi_s\left(\mathbf{s}-\alpha\eta_s\mathbf{d}_s\right),
\]

\[
\rho_{trial}
=\Pi_\rho\left(\rho-\alpha\eta_\rho d_\rho\right).
\]

This is a descent direction whenever either gradient block is nonzero. A
single backtracking factor keeps the joint trial atomic: shift and angle are
accepted together only when the complete loss decreases.

The common `alpha` homothetically shrinks the two block radii. Retain this
joint step for v1 rather than introducing separate block-coordinate line
searches. Diagnostics must nevertheless report whether repeated rejection is
caused predominantly by the shift or angle block; persistent throttling by one
block is evidence for a future block-coordinate or diagonally preconditioned
variant.

### Step-radius parameters

`eta_s` should reuse the existing development parameter `sgd_eta_shift`, but
because the shift gradient block is normalized, the parameter is not an
ordinary gradient-magnitude learning rate. It is the maximum full-step shift
displacement in pixels before backtracking and projection.

A separate developer-only `sgd_eta_inpl` is recommended for `eta_rho`, because
pixels and rotation bins have different units. It is the maximum full-step
continuous-angle displacement in rotation bins. Proposed initial default:

```text
sgd_eta_inpl=0.25 rotation bins
```

This default requires focused parameter-sensitivity validation before it is
accepted. Retain the public `sgd_eta_*` names for consistency with the
development command, but document them in UI/help text as maximum full-step
block displacements. Reusing one value for both blocks would couple unrelated
physical units and is not the preferred mature-optimizer design.

### Per-step procedure

For each direct joint step:

1. Evaluate `L_current` and `grad_current` with
   `gen_raw_euclid_grad_at_angle`.
2. Construct the independently scaled descent proposal and project it into the
   feasible product set

   \[
   C=[-trs,trs]^2\times[r_0-2,r_0+2].
   \]

   With backtracking factor `alpha`, define the actual feasible step

   \[
   \mathbf{p}_\alpha=
   \Pi_C(\mathbf{z}-\alpha\mathbf{D}(\mathbf{g}))-\mathbf{z},
   \]

   where `D(g)` contains the normalized shift block scaled by `eta_s` and the
   signed angle block scaled by `eta_rho`.
3. Stop validly when the projected feasible step is negligible. Do not use the
   raw gradient norm as the constrained-stationarity test: a nonzero gradient
   pointing outside an active shift or angle bound can describe a valid
   constrained solution.
4. Compute the predicted feasible descent `dot_product(grad,p_alpha)`. If it is
   not meaningfully negative at floating-point scale, stop with no numerically
   usable feasible descent. Merely being at a bound is not a stopping rule: a
   component whose descent direction points back into the feasible region must
   remain free to move.
5. Evaluate the full raw loss at the projected trial.
6. Accept only when the trial is finite and satisfies projected Armijo
   sufficient decrease:

   \[
   L(\mathbf{z}+\mathbf{p}_\alpha)
   \leq L(\mathbf{z})+c_1\mathbf{g}^T\mathbf{p}_\alpha,
   \qquad 0<c_1\ll1,
   \]

   Do not impose an independent fixed absolute loss-decrease threshold. Near a
   valid solution, mathematically legitimate improvements can be small; use
   floating-point scale to classify the predicted feasible descent instead.
7. If a speculative trial is non-finite or fails Armijo, reject only that trial,
   halve the common backtracking factor, and retry up to the bounded
   backtracking limit.
8. Stop when no backtracked feasible trial provides sufficient decrease.

The current accepted state is transactional. A non-finite initial/current
accepted evaluation invalidates the optimizer. A non-finite speculative trial
does not: after backtracking is exhausted, return the last accepted finite
state with its existing validity/improvement status.

The optimizer reports a valid evaluation whenever the seed loss and coordinates
are finite and in bounds. It reports improvement only when at least one step
was accepted and the final loss is strictly below the initial loss. These are
separate status values.

Termination diagnostics must distinguish:

- unconstrained stationarity;
- shift-bound-limited stationarity;
- angle-window-limited stationarity; and
- backtracking exhaustion.

An angle-window-limited result with an outward derivative is a valid constrained
result, but not evidence that the underlying continuous pose is stationary.

## Return and rejection contract

The direct joint routine returns the streaming representation:

```text
cxy(1) = -final_raw_loss
cxy(2:3) = final shift in the requested frame
irot = nearest periodic integer rotation bin
rotind_frac = final unwrapped working rotation-grid coordinate
```

`rho_work=rotind_frac` remains unwrapped while optimizing inside
`[r0-2,r0+2]`, which avoids a discontinuity when the local window crosses
rotation index `1` or `nrots`. It must not be written directly as durable
metadata. Finalization produces:

\[
irot=\operatorname{modulo}(\operatorname{nint}(\rho_{work})-1,nrots)+1
\]

and uses the established `store_continuous_e3` convention to persist a
canonical angle in `[0,360)`:

\[
e3=\operatorname{modulo}
\left(360-(\rho_{work}-1)\,dang,360\right).
\]

Optional outputs should mirror `minimize_direct`:

```fortran
accepted_steps
objective_initial
objective_final
evaluation_valid
improved
```

On a numerically invalid initial or current accepted evaluation:

```text
irot = 0
evaluation_valid = false
improved = false
accepted_steps = 0
the strategy preserves the incoming class, rotation, shift, score, and e3
```

On a finite valid seed with no improving trial:

```text
irot = the supplied discrete seed
rotind_frac = real(irot,dp)
cxy(1) = -initial_raw_loss
cxy(2:3) = xy_in in the requested frame
evaluation_valid = true
improved = false
accepted_steps = 0
```

Thus no-improvement is not represented as numerical failure. The caller can
retain the exact selected candidate while still distinguishing a valid flat or
locally stationary objective from an invalid evaluation. This matches the
current conventional joint optimizer's separate validity/improvement semantics
while preserving the stream's `-L` representation.

If one or more steps were already accepted and a later speculative trial is
non-finite, the routine backtracks from the last accepted finite state. If no
finite Armijo trial remains, it returns that last accepted state with
`evaluation_valid=true`, `improved=true`, and the actual positive
`accepted_steps`; it must not roll back to the original seed.

No `exp(-L)` conversion occurs inside `minimize_joint_raw`. The stream keeps
`-L` throughout candidate comparison and converts to SIMPLE's persisted
Euclidean score only at the existing project boundary.

## Strategy integration

The streaming strategy remains a mixed discrete/continuous optimizer:

```text
for each sampled particle
    enumerate classes and rotation bins with gen_raw_euclid_vals
    choose one hard class/rotation winner

    if inpl_cont=no
        refine shift with minimize_direct
    else
        pass that exact winner to minimize_joint_raw
        refine shift and fractional rotation locally without another angle scan

    persist one hard class, nearest rotation bin, continuous e3, and shift
```

Mathematically:

\[
(c^*,r_0^*)=\arg\min_{c,r}L(c,r,\mathbf{s}_0),
\]

followed by local continuous refinement:

\[
(\mathbf{s}^*,\rho^*)
=\operatorname*{local\ argmin}_{\mathbf{s},\rho}
L(c^*,\rho,\mathbf{s}).
\]

This does not make the class variable continuous and does not globally refine
every rotation candidate. It continuously polishes the selected discrete
rotation seed.

### Required strategy state rules

- Dispatch from the existing iteration-local `l_sgd_streaming_active` policy;
  do not infer activity merely from the development command name.
- Construct the existing direct shift objects when SGD is active and
  `inpl_cont=no`.
- Construct the new direct joint object when SGD is active and
  `inpl_cont=yes`.
- Never construct the L-BFGS-B joint object for an SGD-active combined route.
- Continue constructing the existing L-BFGS-B joint object for conventional
  `inpl_cont=yes`, including iterations of `abinitio2D_sgd` where the stage or
  boundary policy has made streaming SGD inactive.
- On direct-joint success, set `best_corr=-L`, update `best_rot`, update the
  shift, and mark `continuous_e3` valid.
- On a valid but non-improving direct-joint evaluation, retain the discrete
  winner, native shift, and `-L` seed merit; reconstruct grid `e3` from the
  integer `inpl`.
- On an invalid direct-joint evaluation, retain the complete incoming selected
  candidate without converting or fabricating a score.
- Store `inpl` as the nearest valid periodic bin and `e3` as the accepted
  real-valued angle. These two fields must describe the same solution.
- Preserve all shift-frame rotations already required by `assign_ori`.

## Commander and UI activation

The commander prohibition is removed only after the numerical routine and
route-identity tests pass.

The final command normalization should be:

```fortran
if( .not. cline%defined('objfun') ) call cline%set('objfun', 'euclid')
if( .not. cline%defined('inpl_cont') ) call cline%set('inpl_cont', 'no')
```

with no rejection of explicit `inpl_cont=yes`. Existing parameter validation
continues to reject continuous refinement for non-Euclidean objectives.

If `sgd_eta_inpl` is approved, expose it only on the developer-visible
`abinitio2D_sgd` descriptor. Do not add it to standard `abinitio2D`.

When activation is implemented, update `doc/policies/abinitio2D_policy.md` in
the same focused change: its present policy correctly says that
`abinitio2D_sgd` rejects `inpl_cont=yes`, which remains true until Phase 4.

## Diagnostics

All new verbose diagnostics remain behind the existing `sgd_diagnostic`
switch. Aggregate counters may be retained without per-particle output.

Recommended direct-joint diagnostics:

- attempted particles;
- accepted particles;
- accepted joint steps;
- rejected/non-finite trials;
- initial and final raw loss summaries;
- absolute fractional-angle increment summaries;
- shift increment summaries;
- bound hits for shift and rotation;
- termination counts for unconstrained stationary, shift-bound-limited,
  angle-window-limited, and backtracking-exhausted states;
- Armijo rejection counts split by non-finite trial, insufficient decrease,
  shift projection, and angle projection;
- best-versus-second-best hard-candidate margin summaries
  `Delta_L=L_second-L_best`, collected with constant memory and no additional
  objective evaluations;
- counts of candidates whose winner margin is within the numerical tie
  tolerance;
- counts showing whether the shift or angle block repeatedly forced the common
  backtracking factor to shrink;
- objective and gradient evaluation counts;
- confirmation that no L-BFGS-B object was used during an SGD-active direct
  joint invocation.

Diagnostic code must not perform vector loss round trips or additional PFTC
evaluations when `sgd_diagnostic=no`. Winner-margin diagnostics retain only the
best and second-best scalar states while the existing candidate stream is
already being traversed; they do not create a top-K table and do not refine the
runner-up in v1.

The winner margin measures the principal scientific approximation in this
design: hard discrete assignment followed by joint continuous pose refinement
of only the selected candidate is cheaper than jointly refining every
class/rotation candidate, but a close runner-up could have a better continuous
minimum. Use the margin distribution during matched workflow validation before
considering a future constant-memory ambiguous-particle extension that refines
both the best and second-best candidates. That extension is explicitly outside
v1 and is not a table-free `prob_snhc` implementation.

## Implementation phases

### Phase 0: freeze the existing behavior

1. Record the reviewed baseline commit and rerun the corrected continuous
   gradient, integer-angle parity, route-identity, SGD base-suite, UI, and
   dispatch tests.
2. Capture the current L-BFGS-B contract: all-angle seed selection at supplied
   shift, `+/-2` local bounds, separate validity/improvement outputs, seed
   return on valid non-improvement, and `irot=0` only for invalid evaluation.
3. Add focused raw-Euclidean seed-equivalence coverage at multiple nonzero
   shifts. Require exact index agreement only for a resolved winner margin;
   near ties require loss equivalence and the documented ordering rule:

   ```text
   argmin(gen_raw_euclid_vals) == select_best_discrete_angle
   ```

4. Confirm the working tree contains no unrelated changes in files selected
   for this work.

Exit criterion: the standalone continuous route has a reproducible baseline.

### Phase 1: private joint-helper refactor

1. Prefer leaving the corrected L-BFGS-B implementation unchanged.
2. Extract only objective-evaluation, projection, periodic-index conversion,
   or shift-frame helpers whose absence would cause real numerical duplication.
3. Do not share seed initialization: conventional joint search owns its
   all-angle rescan; direct joint SGD owns the supplied hard winner.
4. If conventional code changes, keep its loop, tolerances, seed selection,
   return score, and validity/improvement behavior unchanged and rerun the
   focused continuous tests.

Exit criterion: no duplicated PFTC formula and, if conventional code changed,
numerical equivalence within tolerance with no route change.

### Phase 2: isolated direct joint optimizer

1. Add `SHSRCH_JOINT_RAW`, `new_joint_raw`, and
   `grad_shsrch_minimize_joint_raw`.
2. Add a pure projected block-scaled trial/gradient-mapping helper.
3. Implement independently scaled shift/rotation radii, projection, projected
   Armijo backtracking, finite checks, transactional trial rejection, and
   feasible-stationarity termination.
4. Accept the stream's integer winner and native shift explicitly; do not call
   `select_best_discrete_angle`.
5. Return separate `evaluation_valid` and `improved` states.
6. Keep the new route unreachable from production commanders.

Exit criterion: focused numerical tests pass without constructing L-BFGS-B.

### Phase 3: production-context strategy integration

1. Split continuous activation into an SGD-inactive L-BFGS-B route and an
   SGD-active direct-joint route, both restricted to supported raw Euclidean
   mode.
2. Construct the direct joint optimizer instead of the shift-only direct
   optimizer for the combined route; retain shift-only construction for
   `inpl_cont=no`.
3. Dispatch shift-only or direct-joint refinement from the same selected hard
   candidate.
4. Preserve `-L` bookkeeping and orientation metadata.
5. Keep the commander prohibition in place while testing the internal route.

Exit criterion: a production-constructor integration test proves the intended
route and fallback behavior.

### Phase 4: development-command activation

1. Remove the `inpl_cont=yes` rejection from
   `simple_commanders_abinitio2d_sgd.f90`.
2. Update the dispatch test to preserve explicit `inpl_cont=yes`.
3. Add `sgd_eta_inpl` only if approved and describe both `sgd_eta_*` controls
   as maximum full-step block displacements rather than ordinary learning
   rates.
4. Update the development descriptor and `abinitio2D` policy documentation.
5. Do not change standard `abinitio2D` metadata or defaults.

Exit criterion: UI/dispatch tests prove default-off behavior and explicit
combined activation.

### Phase 5: focused and workflow validation

Run the validation matrix below. Do not commit the final activation until the
combined route passes the focused suite and an Oracle/Linux workflow.

## Local Windows/MSYS implementation workflow

The Windows build is a user-local development aid, not an officially supported
SIMPLE platform and not a substitute for the final Oracle/Linux workflow.
Use these repository notes/scripts as the canonical local references:

- `.codex/windows-msys-install-workflow.md`
- `.codex/compile_windows.sh`

Do not launch a bare Windows Bash process. Invoke MSYS2 through its UCRT64
environment wrapper so the compiler, CMake generator, Perl, and Unix tools are
resolved consistently. Use an explicit source and build directory because the
compile script's generic default build path is not named for the `hael_SIMPLE`
checkout:

```powershell
& "C:\msys64\usr\bin\env.exe" `
  MSYSTEM=UCRT64 `
  CHERE_INVOKING=1 `
  /usr/bin/bash -lc '
    set -euo pipefail
    export PATH=/ucrt64/bin:/usr/bin:$PATH
    cd /home/hossainm7/hael_SIMPLE
    SIMPLE_SOURCE_DIR=/home/hossainm7/hael_SIMPLE \
    SIMPLE_BUILD_DIR=/home/hossainm7/hael_SIMPLE_build_debug \
    SIMPLE_BUILD_TYPE=Debug \
    SIMPLE_BUILD_JOBS=12 \
    SIMPLE_CLEAN_BUILD=no \
    SIMPLE_RUN_TESTS=no \
      .codex/compile_windows.sh
  '
```

For the phase-by-phase development loop:

1. Keep `SIMPLE_CLEAN_BUILD=no` for safe incremental compilation. The script
   removes the build tree only when `SIMPLE_CLEAN_BUILD=yes` is explicitly
   supplied, and it refuses unsafe source/root build-directory values.
2. Use at most 12 parallel jobs on the local laptop.
3. Prefer targeted builds after each phase when a full install is unnecessary,
   for example:

   ```bash
   cmake --build /home/hossainm7/hael_SIMPLE_build_debug \
     --target simple_exec simple_test_ui_visibility \
              simple_test_abinitio2d_sgd_dispatch \
              simple_test_sgd_base_suite \
     --parallel 12
   ```

4. Run the numerical continuous-angle tests with their established fixtures
   whenever Phase 1 or Phase 2 touches objective, gradient, projection, or
   angle-conversion code.
5. Inspect `git status --short --branch` after building. Some build variants
   stamp executable sources; generated-only changes must not enter a source
   commit.
6. Use `SIMPLE_CLEAN_BUILD=yes` for a deliberate milestone rebuild before the
   local handoff, not for every edit/compile cycle.
7. Treat successful Windows compilation and focused tests as API, dispatch,
   and numerical-development evidence only. The Oracle/Linux focused suite and
   matched `abinitio2D_sgd` workflow remain required for acceptance.

Known Windows-only test-runner/listing caveats documented in the workflow note
must not be mistaken for scientific validation or silently reported as passed.

### Next-session starting point

Do not begin this direct-solver plan in the next implementation session. Follow
`abinitio2d_sgd_seeded_lbfgsb_plan.md` first and preserve all unrelated
worktree changes. Resume Phase 0 here only after the seeded-L-BFGS-B A/B gate
has measured final loss, evaluation count, and runtime from the same hard seed,
and only after an explicit instruction to continue with solver replacement.

## Validation matrix

### Pure numerical contract

1. Preserve the existing corrected-gradient checks:
   - central differences in `sx`, `sy`, and `rho` at optimized and deliberately
     nonstationary off-grid states;
   - integer-angle loss and `(sx,sy)` gradient parity with
     `gen_raw_euclid_grad_for_rot_8` at several nonzero shifts; and
   - finite results across low/full frequency bands and periodic boundaries.
2. Verify the hard-stream vector winner and
   `select_best_discrete_angle` choose the same integer rotation at the same
   class, particle, and nonzero shift when the winner margin exceeds the
   numerical tie tolerance; for ties/near ties, verify selected-loss
   equivalence and the documented deterministic ordering rule.
3. Verify `minimize_joint_raw` begins from the supplied winner and does not add
   an `nrots`-sized objective scan to its profile counters.
4. Verify projected joint trials never exceed any shift or unwrapped-angle
   bound, including a `+/-2` interval crossing the periodic boundary.
5. Verify projected stationarity at interior, shift-bound, and angle-window
   constrained solutions, including nonzero raw gradients pointing outside
   active bounds.
6. Verify every accepted step is finite, strictly lowers raw loss, and satisfies
   projected Armijo sufficient decrease plus the floating-point guard.
7. Verify an invalid initial/current evaluation preserves all incoming strategy
   state.
8. Verify a non-finite speculative trial triggers backtracking and cannot erase
   earlier accepted finite steps.
9. Verify valid no-improvement returns the supplied seed, `-L_seed`, and
   `accepted_steps=0` without setting `irot=0`.
10. Verify zero gradients and backtracking exhaustion terminate safely.
11. Verify unwrapped accepted `rho` maps to the expected periodic integer bin
    and canonical `[0,360)` `e3`, including both periodic boundaries.
12. Verify the return value is `-L`, not `exp(-L)`.

### Orthogonality and route identity

Prove every route state in the compatibility table:

- conventional default does not construct either joint optimizer;
- conventional `inpl_cont=yes` constructs only the established L-BFGS-B
  joint optimizer;
- SGD default constructs only the direct shift optimizer;
- SGD plus `inpl_cont=yes` on an SGD-inactive iteration constructs the
  established L-BFGS-B joint optimizer;
- SGD plus `inpl_cont=yes` on an SGD-active iteration constructs the direct
  joint optimizer and never constructs or calls L-BFGS-B.

Exercise at least one transition between an inactive/warmup iteration and an
active streaming iteration. This proves that continuous refinement changes
optimizer with the existing SGD schedule rather than changing the schedule
itself.

Also verify that the conventional route still performs its owned all-angle
seed selection, whereas the direct joint route consumes the already-selected
stream seed without a second scan.

The test must use production constructors, not test-only state injection.

### Synthetic truth recovery

Extend the deterministic V4-style test over several:

- classes;
- discrete rotation seeds;
- off-grid true angles;
- positive and negative shifts;
- noise levels;
- periodic angle boundaries.

Report:

- class accuracy;
- nearest rotation-bin accuracy;
- continuous angle error;
- shift residual;
- objective reduction;
- accepted-step count;
- bound compliance;
- no-improvement preservation;
- projected-stationarity termination reason; and
- hard-winner margin `L_second-L_best`.

### Regression suites

At minimum:

```text
simple_test_ui_visibility
simple_test_abinitio2d_sgd_dispatch
simple_test_continuous_inplane_rotation2D_route_identity
simple_test_continuous_inplane_rotation2D_stage1_validation
simple_test_sgd_base_suite
```

The standalone continuous tests must pass without weakened tolerances. Extend
the route-identity or SGD suite rather than creating a test-only numerical
implementation of the joint loss.

### Matched workflow comparison

Fork all modes from the same preprocessing/checkpoint state:

1. conventional baseline;
2. conventional continuous in-plane;
3. shift-only streaming SGD;
4. direct-joint streaming SGD.

Compare:

- normal completion;
- runtime;
- best and median resolution;
- classes below the chosen resolution threshold;
- ranked and active classes;
- zero-support classes;
- class-population distribution;
- angle/shift increment summaries;
- objective reduction;
- hard-winner margin distribution and near-tie fraction;
- frequency of angle-window-limited solutions;
- evidence that the common backtracking factor is or is not repeatedly
  throttled by one coordinate block;
- proof of optimizer route.

Use more than one dataset before selecting a default or recommending the
combined mode scientifically. `inpl_cont` remains default-off regardless of
the result of this development comparison.

## Commit boundaries

If implementation is approved, keep commits reviewable:

1. private joint-helper refactor plus equivalence test;
2. isolated direct joint optimizer plus numerical/seed-equivalence tests;
3. strategy integration plus production-constructor route-identity test;
4. development-command/UI/policy activation plus dispatch test;
5. workflow scripts or scientific reporting, if changed.

Do not mix continuous-refinement feature development, standard `abinitio2D`
changes, or unrelated staged work into these commits.

## Approved design decisions

The following decisions were approved on 2026-08-06:

1. Use the type-bound name `minimize_joint_raw` and implementation name
   `grad_shsrch_minimize_joint_raw`.

   `minimize_direct_joint` would describe the direct first-order algorithm more
   explicitly, but `minimize_joint_raw` correctly emphasizes the API
   contract that matters to its stream caller: it consumes and returns the
   raw Euclidean loss representation rather than the established
   `exp(-L)` joint score. The dedicated `SHSRCH_JOINT_RAW` mode, constructor,
   comments, and route-identity tests must make the absence of L-BFGS-B
   unambiguous.

2. Add a separate developer-only `sgd_eta_inpl` parameter.

   `sgd_eta_shift` remains active and is not made obsolete. It continues to
   control pixel-space steps in the existing shift-only stream and controls
   the `(sx,sy)` block in the direct joint stream. `sgd_eta_inpl` controls the
   fractional-rotation-bin step independently. Under normalized block descent,
   these are maximum full-step trust radii/displacements, not conventional
   gradient-magnitude learning rates. Unused SGD variables may be audited only
   after the combined route is implemented and validated.

3. Retain the established local angle window of `+/-2` rotation bins.

   This bounds the continuous polish around the selected discrete seed and is
   substantially cheaper than a global continuous-angle search. The discrete
   selector has already scanned the rotation grid, so the joint optimizer is
   responsible only for local off-grid refinement.

4. Refine only the selected discrete class/rotation candidate.

   This requires one joint optimization per sampled particle rather than one
   optimization for every class/rotation candidate. It preserves the intended
   table-free hard-assignment cost. The accepted tradeoff is that a different
   discrete candidate might have won if every candidate had received an
   expensive continuous refinement.

These decisions complete the long-term design gates. Source implementation is
deferred until the seeded-L-BFGS-B current plan has been implemented, tested,
and reviewed; resuming this plan still requires an explicit instruction.

## Repository-driven plan update (2026-08-07)

The continuous-gradient correction in `1f2b81f40` resolves the numerical
prerequisite that originally blocked this design. The plan now treats
`gen_raw_euclid_grad_at_angle` as ready for reuse, subject to rerunning its
existing tests in the implementation branch. It does not propose another
gradient formula.

The same repository update also clarified that conventional and streaming
joint search must not be made identical at the orchestration level:

- conventional `minimize_joint` deliberately reselects the best grid angle at
  the supplied shift and uses L-BFGS-B;
- proposed `minimize_joint_raw` deliberately accepts the stream's existing
  hard winner and uses bounded block-scaled first-order descent with projected
  Armijo backtracking; and
- a focused seed-equivalence test connects those initial states under the raw
  Euclidean objective without paying for a redundant production scan.

The replacement is therefore iteration-local: `abinitio2D_sgd` can retain
established L-BFGS-B refinement during its conventional warmup/inactive
iterations, then switch to `minimize_joint_raw` exactly when the existing stage
policy activates the hard stream. The new optimizer does not alter that stage
policy; it supplies the continuous optimizer used by its active SGD branch.

Finally, the direct routine will expose validity and improvement separately.
Only invalid evaluation uses `irot=0`; valid no-improvement returns the exact
seed in raw-loss representation. These updates supersede the earlier draft's
assumption that every non-improving joint result should be treated as failure.

## Theoretical-review refinements (2026-08-07)

The independent plan review was incorporated without changing the approved
discrete/continuous architecture:

1. The inner three-variable method is defined as bounded block-scaled
   projected-gradient descent with Armijo backtracking. Its step has a
   first-order block trust-region interpretation, while the outer
   sampled-particle workflow remains stochastic.
2. Convergence uses a projected feasible-step/gradient-mapping criterion rather
   than the raw gradient norm.
3. Backtracking first requires meaningfully negative predicted feasible
   descent and then uses projected Armijo sufficient decrease; it does not add
   a fixed absolute loss-decrease threshold.
4. The redundant `rotind_seed` argument is removed from v1; the hard integer
   winner initializes `rho=real(irot,dp)`.
5. A non-finite speculative trial is transactional and cannot invalidate or
   erase the last accepted finite state.
6. Seed-equivalence testing is winner-margin and tie aware.
7. Optimization keeps `rho` unwrapped locally, while durable `e3` is
   canonicalized through the established `[0,360)` storage convention.
8. Scientific validation records the best/second-best hard-candidate margin
   with constant memory. Refining a runner-up, adding a top-2 ambiguity policy,
   or implementing a table-free `prob_snhc` analogue remains future work and
   is not part of v1.

The v1 joint step remains atomic with one common backtracking factor. A
block-coordinate or diagonally/spectrally preconditioned variant is considered
only if diagnostics show that one block consistently throttles progress or the
simple method converges too slowly.
