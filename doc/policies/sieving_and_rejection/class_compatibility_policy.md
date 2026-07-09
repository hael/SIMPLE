# Class Compatibility Policy

This document defines the behavioral policy for class-average size
compatibility in SIMPLE, implemented by the module
[../../src/main/class/simple_class_compatibility.f90](../../src/main/class/simple_class_compatibility.f90).

The policy governs how support-model fitting, inference, convergence,
and telemetry behave. It is intended to make rejection decisions stable,
explainable, and safe to extend.

## 1. Scope

The class compatibility workflow is responsible for:

1. preprocessing class averages into binary support masks;
2. extracting Feret min/max dimensions per class average;
3. fitting latent support bounds `(c, b, a)` from selected classes;
4. rejecting size-incompatible classes during inference;
5. exposing fit metrics and convergence state for monitoring.

It is not responsible for quality-feature extraction, class-average quality
modeling, or queue/chunk orchestration. Those remain owned by sieve and
quality modules.

## 2. Public Contract

Public surface (type-bound methods on `class_compatibility`):

1. `new`, `kill`, `kill_support_model`
2. `train` via generic (`train_1` on project, `train_2` on stack)
3. `infer`
4. `converged`
5. `get_support_model_metrics`

Public metrics type:

- `support_model_metrics` with fields:
  - `axis_c`, `axis_b`, `axis_a`
  - `delta_c`, `delta_b`, `delta_a`
  - `delta_valid`, `valid`, `converged`

## 3. Core Semantics

### 3.1 Validity

`valid` means only one thing: the module currently has a fitted axis model.

- `valid=.false.` when no successful fit has occurred, or after
  `kill_support_model`.
- `valid` must not be used as a proxy for training batch presence.

### 3.2 Convergence

Convergence is axis-stability based and requires a previous fit.

Given previous axes `(c_prev, b_prev, a_prev)` and current axes `(c, b, a)`,
axis deltas are:

- `dc = abs(c - c_prev)`
- `db = abs(b - b_prev)`
- `da = abs(a - a_prev)`

An axis is stable when either absolute or relative tolerance passes:

- `d <= CONV_ABS_EPS`
- or `d <= CONV_REL_EPS * max(abs(prev), abs(curr), 1.0)`

`converged=.true.` only when:

1. a previous valid fit exists, and
2. all three axes are stable.

### 3.3 Delta Validity

`delta_valid` indicates that delta values compare two real fits.

- First successful fit: `delta_valid=.false.`.
- Subsequent successful fits: `delta_valid=.true.`.

## 4. Training Policy

Training data is appended as batches (`training_set(:)`), each batch storing
vector pairs (`min_dim(:)`, `max_dim(:)`) extracted from selected classes.

Policy requirements:

1. each appended batch must satisfy `size(min_dim) == size(max_dim)`;
2. batches with fewer than 3 selected classes are ignored;
3. update must preserve prior valid fit if no new successful fit is found;
4. batch weighting must be balanced so each batch contributes equal total
   weight regardless of batch size.

Fitting strategy policy:

1. grid search over `relax`, `qlo`, `qhi` candidates;
2. stage 1: choose `b` maximizing weighted in-interval support;
3. stage 2: estimate `c` and `a` from weighted quantiles on the support set;
4. enforce axis ordering invariant `c <= b <= a` via clamped fallback;
5. score preference is: support first, compactness second, lower relaxation
   third.

## 5. Inference Policy

Inference is state-preserving unless a class is proven incompatible.

For each active class:

1. preprocess image and compute `(min_dim, max_dim)`;
2. if model is invalid, skip compatibility rejection;
3. test compatibility against relaxed `b` interval and relaxed global `c/a`
   bounds;
4. apply edge-rescue slack when strict compatibility fails;
5. incompatible class must be set to `state=0` and receive rejection reason
   `class_compatibility: size_incompatible_subset`.

After per-class decisions, selection must be propagated through
`map_cavgs_selection`.

## 6. Preprocessing Policy

Preprocessing must remain deterministic and minimal:

1. edge-average removal + low-pass bandpass
2. resize to fixed working box (`PREPROCESS_BOXSIZE`)
3. Otsu segmentation
4. fixed morphological closing depth (`PREPROCESS_MORPH_SIZE`)
5. largest connected component extraction
6. Feret min/max measurement on final mask

Changes to preprocessing parameters are policy-significant and require
re-baselining tests and downstream acceptance criteria.

## 7. Observability Requirements

Callers (notably sieve orchestration) should log:

1. current axes `a/b/c`
2. deltas `da/db/dc`
3. `valid`, `delta_valid`, `converged`

Convergence transitions should be logged explicitly when first observed.

## 8. Failure Handling

The module must fail fast on structural inconsistencies:

1. image/state size mismatch
2. malformed training batch allocation state
3. min/max vector length mismatch

Non-fatal conditions (no fit update) must return cleanly without corrupting
existing valid model state.

## 9. Test Policy

Policy-level tests must cover:

1. default state and lifecycle reset behavior
2. first-fit validity semantics (`valid=true`, `delta_valid=false`)
3. repeat-fit convergence behavior (`delta_valid=true`, stable deltas)
4. infer behavior with invalid model (no forced rejection)
5. infer behavior with trained model (outlier rejection)
6. small-batch no-fit path
7. kill/reset semantics after valid fit

Reference tester module:
[../../src/main/class/simple_class_compatibility_tester.f90](../../src/main/class/simple_class_compatibility_tester.f90).

## 10. Change Checklist

When modifying class compatibility behavior:

1. preserve `valid` as fitted-axis validity only;
2. preserve axis ordering invariant `c <= b <= a`;
3. preserve convergence gate requiring prior fit and 3-axis stability;
4. keep rejection reason string stable unless migration is intentional;
5. update tester coverage for any policy-level behavior change;
6. update this policy document in the same change when contract changes.
