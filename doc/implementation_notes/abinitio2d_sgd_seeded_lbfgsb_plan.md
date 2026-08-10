# Hard-stream-seeded L-BFGS-B for `abinitio2D_sgd`

## Status

Current incremental implementation plan as of 2026-08-10.

Repository baseline reviewed: `master` at `c769408f1` (`more
inpl_cont=yes updates`). This plan must be implemented and evaluated before the
longer-term direct-solver direction in
`abinitio2d_sgd_direct_joint_optimizer_plan.md` is resumed.

This note authorizes no default behavior change. Standard `abinitio2D` and the
default `abinitio2D_sgd inpl_cont=no` route remain unchanged throughout the
experiment. Remove the current `abinitio2D_sgd inpl_cont=yes` commander guard
only after the isolated seeded route and production-context tests pass.

## Objective

Give SIMPLE's established three-variable L-BFGS-B joint optimizer the class,
integer rotation, and shift seed already selected by the table-free hard
stream:

```text
hard-stream winner -> seeded L-BFGS-B over (sx,sy,rho)
```

The immediate question is whether eliminating the wrapper-level all-angle
rescan provides a useful and scientifically safe improvement while retaining
the mature continuous solver. This experiment deliberately separates:

1. the value of reusing the stream-owned discrete seed; from
2. the value, if any, of replacing L-BFGS-B itself.

The future projected direct optimizer will consume the same seed/problem
contract. Hard-stream seed reuse therefore benefits both solvers; it is not a
feature specific to projected descent.

## Important terminology

The existing all-angle rescan is not part of Algorithm 778/L-BFGS-B. It is
performed by SIMPLE's `grad_shsrch_minimize_joint` wrapper before the generic
optimizer starts:

```text
current minimize_joint
    select_best_discrete_angle(xy_in)
    set rho bounds to selected_irot +/- 2
    initialize (xy_in,selected_irot)
    call L-BFGS-B
    validate and finalize
```

This plan does not edit or progressively replace the generalized Cauchy point,
free-variable subspace solve, BFGS update, or safeguarded Wolfe line search in
`simple_opt_lbfgsb.f90`. Those components are coupled parts of a mature solver.
The plan changes seed ownership at the PFTC wrapper boundary.

## Current route facts

The table-free active stream already:

1. evaluates raw Euclidean losses for every eligible class/rotation at one
   current native-frame shift;
2. selects one hard class/rotation winner with lower raw loss preferred;
3. stores the internal merit as `-L` so existing larger-is-better strategy
   bookkeeping remains valid; and
4. runs fixed-angle `minimize_direct` for the selected winner.

The conventional `inpl_cont=yes` route fixes the selected class but discards
the incoming rotation, rescans every rotation at `xy_in`, then invokes
L-BFGS-B over:

\[
z=(s_x,s_y,\rho),
\]

with shift limits and an unwrapped local angular interval of the newly selected
rotation plus or minus two bins. Its objective and gradient are supplied by
`gen_raw_euclid_grad_at_angle`.

The hard stream and conventional selector are expected to choose the same
rotation for the same class, particle, and native shift because
`exp(-L)` is monotone. Floating-point ties and near ties require a margin-aware
test rather than unconditional bitwise index equality.

## Compatibility contract

| Command and iteration state | Required behavior |
| --- | --- |
| `abinitio2D`, `inpl_cont=no` | Unchanged |
| `abinitio2D`, `inpl_cont=yes` | Existing rescan-seeded L-BFGS-B, unchanged |
| `abinitio2D_sgd`, `inpl_cont=no` | Existing hard stream plus fixed-angle `minimize_direct`, unchanged |
| `abinitio2D_sgd`, `inpl_cont=yes`, SGD inactive/warmup | Existing conventional route selected by current stage policy |
| `abinitio2D_sgd`, `inpl_cont=yes`, SGD active | Hard-stream winner plus seeded L-BFGS-B, with no second rotation scan |
| Terminal conventional/probability pass | Existing route and artifact contract |

Additional invariants:

- `inpl_cont=no` remains the default.
- No new public parameter is needed for this experiment.
- No SoftMax, top-K, `eulprob_tab2D`, or likelihood-table machinery is added to
  the active hard stream.
- Stage activation, mini-batch sampling, zero-support recovery, and terminal
  pass policy are unchanged.
- Class selection remains fixed during the local continuous solve.
- The local angular window remains `irot_seed +/- 2` bins.
- A valid non-improving solve preserves the supplied seed exactly.
- A numerically invalid solve preserves the complete incoming strategy state.
- Standard `abinitio2D` must not construct or call the seeded entry point.

## Shared joint-pose boundary

Do not begin by creating a broad new public module. The existing numerical
owner, `simple_pftc_shsrch_grad`, is already a Fortran module and owns the
builder/PFTC context, objective callbacks, optimizer object, limits, frame
conversion, and profiling counters.

Introduce the smallest private seed/problem helpers needed to prevent the two
entry points from drifting. Candidate responsibilities are:

```fortran
prepare_joint_seed(irot_seed, xy_seed, joint_lims, x0)
    validate integer rotation and finite native shift
    set x0 = [xy_seed,real(irot_seed)]
    set rho bounds = real(irot_seed) +/- 2

evaluate_joint_raw(...)
    call gen_raw_euclid_grad_at_angle
    maintain objective/gradient counters

validate_joint_solution(...)
    distinguish finite valid evaluation from strict improvement
    verify all returned coordinates are inside bounds

finalize_joint_solution(...)
    preserve the exact seed on valid no-improvement
    map improved rho to the nearest periodic integer bin
    retain unwrapped fractional rho for canonical e3 storage
    convert shift frame only when requested
```

The conventional initializer remains:

```text
select_best_discrete_angle -> prepare_joint_seed -> L-BFGS-B
```

The new initializer is:

```text
supplied hard winner -> prepare_joint_seed -> L-BFGS-B
```

Prefer leaving the established `minimize_joint` source unchanged during the
first isolated experiment. A small amount of temporary duplication is safer
than refactoring the accepted route before A/B equivalence is known. Extract
shared preparation/finalization helpers only after focused tests prove the new
entry point's contract.

## Proposed seeded API

Add a dedicated type-bound entry point whose name states both seed ownership
and solver identity, for example:

```fortran
procedure :: minimize_joint_seeded_lbfgsb
```

Its logical contract is:

```fortran
function minimize_joint_seeded_lbfgsb(self, irot_seed, xy_seed, sh_rot, &
    irot_out, rotind_frac, evaluation_valid, improved, &
    objective_initial, objective_final) result(cxy)
```

Exact argument order may follow local Fortran style, but these ownership rules
are mandatory:

- `irot_seed` is `intent(in)` and is never overwritten.
- `xy_seed` is the native-frame shift at which the hard stream evaluated the
  winning class/rotation.
- `irot_out` is separate from the seed and is the nearest periodic bin of the
  accepted fractional result.
- `rotind_frac` is optimizer-local/unwrapped until final metadata conversion.
- The routine must not call `select_best_discrete_angle`,
  `gen_raw_euclid_vals`, `gen_objfun_vals`, or any all-angle evaluator.
- The routine must call the existing L-BFGS-B optimizer through the same joint
  raw objective callback and use the same `factr`, `pgtol`, and evaluation
  limit as the conventional route.

For stream compatibility, define the dedicated result representation as:

```text
cxy(1)   = -final_raw_loss
cxy(2:3) = accepted shift in the requested frame
```

Do not compute `log(exp(-L))` to recover the raw merit. The objective callback
already provides finite raw `L`; use it directly and convert to SIMPLE's
established `exp(-L)` project score only at the existing assignment boundary.

The existing `minimize_joint` continues returning its established `exp(-L)`
representation. The two APIs intentionally differ at their caller boundary,
so names, comments, and tests must make the representation explicit.

## Validity and rejection contract

The seeded L-BFGS-B wrapper must preserve the mature solver while making its
transaction boundary explicit.

On an invalid seed or invalid initial objective:

```text
evaluation_valid = false
improved = false
irot_out = 0
strategy preserves incoming class, rotation, shift, score, and e3
```

On a finite valid seed with no strict improvement:

```text
evaluation_valid = true
improved = false
irot_out = irot_seed
rotind_frac = real(irot_seed,dp)
cxy(1) = -objective_initial
cxy(2:3) = xy_seed in the requested frame
```

On improvement:

```text
evaluation_valid = true
improved = true
cxy(1) = -objective_final
irot_out = nearest periodic bin of rotind_frac
```

An L-BFGS-B task that stops because of its evaluation limit is not
automatically invalid. Mirror the existing wrapper: accept a finite in-bounds
result only when it strictly improves the initial raw loss within the existing
floating-point-scaled tolerance; otherwise preserve the seed.

## A-versus-B experiment

The first scientific comparison changes seed ownership only:

| Case | Seed preparation | Continuous solver |
| --- | --- | --- |
| A | `select_best_discrete_angle(xy_seed)` | Existing L-BFGS-B |
| B | Supplied hard-stream `(class,irot_seed,xy_seed)` | Same L-BFGS-B |

For every comparison, use the same:

- PFTC context and frequency limits;
- particle and class reference;
- native-frame shift;
- raw Euclidean objective;
- shift bounds and plus-or-minus-two-bin angle window;
- L-BFGS-B tolerances and evaluation limit.

First calculate the hard winner and its best-versus-second-best margin. When
the margin is resolved beyond scaled numerical tolerance, require A's rescan to
select the same integer bin. For a tie or near tie, require selected raw-loss
equivalence and record the deterministic candidate-order difference rather
than failing on the index alone.

When the selected seed agrees, compare:

- identical initial raw loss;
- final raw loss within numerical tolerance;
- validity and improvement outcome;
- final shift and fractional angle, interpreted with periodicity;
- bound compliance;
- objective and gradient evaluations; and
- wall time over a sufficiently large particle sample.

Do not require bitwise-identical final coordinates in a shallow local valley.
The principal numerical gate is comparable final raw loss from an identical
seed. Profile accounting must prove that B does not add the conventional
`nrots`-sized scan.

## Implementation phases

### Phase 0: freeze current behavior

1. Record the current commit and dirty worktree before editing.
2. Run the focused corrected-gradient, integer-angle parity, route-identity,
   UI, dispatch, and SGD base-suite tests.
3. Add or extend a focused hard-winner/rescan seed-equivalence test at zero and
   nonzero shifts with tie-aware assertions.
4. Capture current L-BFGS-B initial/final loss and profile counters.

Exit criterion: reproducible Case-A baseline and verified hard winner.

### Phase 1: isolated seeded L-BFGS-B entry point

1. Add `minimize_joint_seeded_lbfgsb` in
   `simple_pftc_shsrch_grad.f90`.
2. Reuse the existing joint raw objective callback and optimizer object.
3. Skip every all-angle evaluator.
4. Implement the explicit raw-merit and validity/no-improvement contracts.
5. Keep the method unreachable from production strategies and commanders.

Exit criterion: the new API compiles and focused invalid, no-improvement,
improvement, boundary, and periodic-window tests pass.

### Phase 2: focused A/B validation

1. Run Cases A and B from identical hard seeds across multiple classes,
   rotations, shifts, and boundary-crossing angle windows.
2. Compare raw losses, coordinates, outcomes, and profile counters.
3. Confirm B removes exactly the wrapper-level rotation rescan without changing
   the L-BFGS-B solver configuration.
4. Record median and tail objective/gradient evaluation counts and runtime.

Exit criterion: seeded L-BFGS-B is numerically equivalent from resolved common
seeds and the removed scan is demonstrated by counters, not log absence.

### Phase 3: production-context integration behind the existing guard

1. Split continuous route state into conventional and active-stream ownership;
   do not overload one Boolean with both meanings.
2. During an SGD-active combined route, preserve the hard stream's selected
   class, integer rotation, and exact native shift.
3. Construct the existing joint L-BFGS-B object and call only the seeded entry
   point after hard selection.
4. Retain fixed-angle `minimize_direct` for `inpl_cont=no`.
5. Retain the established conventional route during SGD-inactive/warmup and
   terminal passes.
6. Keep `abinitio2D_sgd inpl_cont=yes` rejected by the commander while internal
   constructor/strategy tests run.

Exit criterion: a production-context test proves route construction, supplied
seed preservation, raw-score semantics, and absence of the second scan.

### Phase 4: development-command activation

1. Remove only the `abinitio2D_sgd inpl_cont=yes` rejection.
2. Update UI/dispatch tests so explicit `inpl_cont=yes` reaches the development
   command unchanged.
3. Add no new public switch and change no standard `abinitio2D` descriptor.
4. Verify stage-policy transitions: inactive iterations use the conventional
   route; active hard-stream iterations use seeded L-BFGS-B.

Exit criterion: default-off and explicit opt-in route tests pass.

### Phase 5: workflow evidence and decision

1. Run the focused suite on Oracle/Linux.
2. Run matched `abinitio2D_sgd` workflows from the same extracted project and,
   where possible, the same pre-SGD checkpoint.
3. Record raw-loss outcomes, scan/evaluation counts, runtime, resolution, class
   occupancy, zero-support recovery, and metadata validity.
4. Decide whether seeded L-BFGS-B becomes the retained active joint route.
5. Only after this decision consider resuming the projected-direct global plan.

Exit criterion: evidence attributes any gain to seed reuse separately from a
future solver replacement.

## Validation matrix

### Numerical/API tests

- invalid integer seed;
- non-finite shift seed;
- initial objective finite/non-finite;
- valid no-improvement exact seed preservation;
- strict improvement;
- shift-bound solution;
- angle-window-bound solution;
- plus-or-minus-two window crossing the periodic grid boundary;
- nearest integer bin and canonical e3 consistency;
- native/reference shift-frame round trip;
- raw `-L` result sign and project-boundary `exp(-L)` conversion;
- objective/gradient profile balance; and
- zero `nrots`-sized scan inside the seeded call.

### Route tests

- standard `abinitio2D` never calls seeded L-BFGS-B;
- standard continuous `abinitio2D` retains its rescan-seeded L-BFGS-B;
- `abinitio2D_sgd inpl_cont=no` retains direct shift-only refinement;
- SGD-inactive combined iterations retain the conventional route;
- SGD-active combined iterations use the supplied hard winner exactly once;
- invalid seeded refinement preserves the incoming hard assignment; and
- terminal conventional behavior remains unchanged.

### Scientific comparison

For Cases A and B report at minimum:

```text
same particle/class/native shift
hard winner and winner margin
rescan winner
initial raw loss
final raw loss
final sx, sy, rho
objective evaluations
gradient evaluations
termination/convergence state
bound status
wall time
```

Aggregate over enough particles to report median and tail behavior. A single
truth-controlled particle proves mechanics but cannot establish runtime or
scientific equivalence.

## Windows/MSYS development workflow

Use the repository-local development references:

- `.codex/windows-msys-install-workflow.md`
- `.codex/compile_windows.sh`

Launch UCRT64 Bash through `C:\msys64\usr\bin\env.exe`, set an explicit
`SIMPLE_BUILD_DIR`, use at most 12 local jobs, and keep incremental builds
non-clean during development. Prefer targeted builds of `simple_exec`, focused
continuous-route tests, UI/dispatch tests, and `simple_test_sgd_base_suite`.

Windows compilation is an early API/numerical gate. Final acceptance requires
the focused Oracle/Linux suite and matched workflow evidence.

## Commit boundaries

Keep implementation commits reviewable:

1. focused seed-equivalence and Case-A baseline tests;
2. isolated seeded-L-BFGS-B numerical API and tests;
3. production-context strategy integration;
4. commander/UI activation and route tests; and
5. workflow/checker updates, only if needed.

Do not mix unrelated 3D, OpenMP-offload, generated, or local workflow changes
into these commits.

## Decision gate for the global direct plan

After seeded L-BFGS-B is measured:

- If it removes the rescan and retains quality at useful lower runtime, it is a
  defensible endpoint even if the direct solver is never implemented.
- If the future projected solver is tested, compare it against Case B from the
  identical hard seed, not against Case A with the extra scan.
- Replace L-BFGS-B only if the projected method reaches comparable final raw
  losses and improves evaluation count, runtime, robustness, or maintainability
  in measured production use.
- If no short- or medium-term benefit is demonstrated, retain seeded L-BFGS-B
  and do not claim solver replacement as an improvement.
