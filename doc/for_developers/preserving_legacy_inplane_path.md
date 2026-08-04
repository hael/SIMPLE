# Preserving the Legacy In-Plane Path

## Purpose

The continuous-angle mathematics was not the critical problem in the original
implementation. The critical problem was the placement of the feature boundary.

`inpl_cont=no` means:

> Run the historical in-plane algorithm without any continuous-angle behavior.

What it actually meant was:

> Do not run the final joint optimizer.

The shared Euclidean shift optimizer had already been changed to use the
continuous angle path, regardless of the value of `inpl_cont`. As a
result, the default-off route was no longer the legacy route.

## Intended invariant

For an opt-in numerical development, the following must hold:

```text
option=off  =>  historical construction, control flow, numerical operations,
                 metadata, and outputs
```

This is stronger than saying that the new final stage is skipped. Every earlier
operation shared with the new development must also retain its historical
configuration.

For this development, `inpl_cont=no` must guarantee all of the
following:

- the historical two-shift L-BFGS-B object was constructed unchanged;
- the historical discrete angle update was used;
- no continuous angular coefficients or Newton updates were evaluated;
- no joint `(sx,sy,theta)` optimizer was constructed;
- shift rotate-back used the selected discrete angle;
- `e3` was reconstructed from the integer `inpl` value; and
- probabilistic and final class-average stages saw only the historical pose.

## Where the implementation crossed the boundary

### 1. The public option gated only the joint stage

At the strategy level, the original implementation derived a flag equivalent
to:

```fortran
use_joint_angle = inpl_cont == 'yes' .and. euclidean_is_active
```

That flag controlled construction and invocation of the three-variable joint
optimizer. This part was gated.

However, the ordinary shift-search objects were constructed identically for
`inpl_cont=yes` and `inpl_cont=no`. The option was not passed into the
lower-level object that owned the angle callback.

The strategy therefore protected only the last new step, not the earlier
behavior that had also changed.

### 2. Objective selection was used as feature selection

Inside `pftc_shsrch_grad`, the angle-update routine selected its algorithm using
the objective type:

```fortran
if( pftc%is_euclid_objfun() )then
    ! continuous angular coefficients and Newton updates
else
    ! discrete angle update
endif
```

The main minimization routine used the same test when initializing the angle,
updating it after shift optimization, interpreting the objective value, and
rotating the final shift.

Consequently, every ordinary Euclidean search took the continuous callback
path. The public option was never consulted at the point where the numerical
algorithm actually changed.

This is the central mistake:

> “Euclidean objective” is a capability condition, not an activation condition.

Euclidean data makes continuous refinement possible. It does not mean the user
requested continuous refinement.

The required condition was conceptually:

```text
use continuous callback = requested callback mode
                          AND raw Euclidean capability
```

not simply:

```text
use continuous callback = Euclidean objective
```

### 3. The default construction silently acquired new behavior

The low-level shift-search object had no explicit property distinguishing the
legacy discrete angle update from the continuous angle callback. Its default
construction therefore no longer represented the legacy algorithm for a
Euclidean objective.

This made it impossible for the high-level `no` route to preserve legacy
behavior merely by skipping the joint object. The shared object had already
changed underneath it.

### 4. The capability check was initially too broad

The continuous derivative was implemented for the raw Euclidean objective, but
the original gateway tested only whether the objective family was Euclidean.
That also admitted the hybrid/denoised Euclidean mode, for which a matching
continuous derivative had not been implemented.

The activation boundary therefore needed both:

```text
objfun=euclid AND .not. l_objfun_den
```

This is a separate issue from the default-off regression, but it has the same
root cause: capability and activation were inferred too indirectly.

## Why the tests did not catch it

The original focused test was strong on numerical identity of the new
coefficient formulation. It compared:

- the legacy Euclidean score;
- the raw FFT loss;
- the scalar loss/gradient path;
- the continuous angular coefficient series; and
- continuous angular derivatives.

Those tests answered:

> Is the new continuous formulation mathematically consistent with the
> Euclidean grid objective?

They did not answer:

> Does `inpl_cont=no` avoid the continuous formulation entirely?

There was no route-identity test that constructed the strategy twice, once with
the option off and once with it on, and inspected which low-level algorithms
were active. There was also no default-off final-metadata test requiring every
`e3` value to remain consistent with its integer `inpl` value.

The missing test was the negative test: proving that the new code does not run
when it has not been requested.

## Corrected design

The production interface is:

```text
inpl_cont=yes|no
```

The route is resolved before the numerical objects are constructed.

### `inpl_cont=no`

- Constructs the historical shift optimizer.
- Uses the historical discrete angle update.
- Does not construct either continuous implementation.
- Does not create or propagate continuous pose state.
- Reconstructs `e3` from the integer rotation index.

### `inpl_cont=yes`

- Explicitly constructs the three-variable object with `new_joint`.
- Applies joint refinement only after discrete candidate selection.

The active route additionally requires the raw Euclidean capability. An
ineligible request resolves to effective mode `no`.

The low-level properties default to false. Thus default construction once again
means legacy behavior; continuous behavior cannot be inferred from the
objective type alone.

## Required regression tests

The route-identity test should assert the construction itself:

| Mode | Joint optimizer | Legacy seed angle update |
|---|---:|---:|
| `no` | off | on |
| `yes` | on | off |

It should also verify:

- a hybrid/denoised request resolves to effective `no`;
- probabilistic stages retain the requested active implementation;
- `no` produces only grid-consistent `e3` values;
- `yes` can produce continuous `e3` values; and
- the `abinitio2D` controller propagates the option through every child and
  final stage.

A fixed-seed, full-workflow comparison remains the strongest confirmation of
legacy identity:

```text
historical executable == current executable with inpl_cont=no
```

The comparison should include assignments, shifts, `inpl`, `e3`, class
averages, and convergence history—not only whether the program completes.

## General lesson

When adding an opt-in numerical algorithm:

1. Define the off-state invariant before implementing the new method.
2. Carry the activation state down to the lowest layer where behavior changes.
3. Keep capability checks separate from activation checks.
4. Make legacy behavior the constructor default.
5. Construct new numerical objects only inside the active route.
6. Keep new metadata invalid unless an active route explicitly produces it.
7. Test the negative path before testing the new path.
8. Verify the final consumer, not only the local optimizer.

The practical review question is:

> If the option is off, is there any line in the new numerical implementation
> that can still execute?

In the original implementation, the answer was yes: the objective-type branch
activated continuous angle refinement inside the shared optimizer. That is why
the legacy path was not preserved.
