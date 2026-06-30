# Abinitio3D SGD Refinement Implementation Plan

## Summary

This note describes how to turn SIMPLE's existing `abinitio3D`
importance-sampled fractional-update machinery into explicit stochastic
gradient descent machinery while preserving the two implementation assets that
matter most:

- Kaiser-Bessel convolution interpolation for Fourier-domain reconstruction and
  reprojection
- fast polar Fourier matching for particle-domain pose, shift, state, and
  candidate assignment

The current `abinitio3D` pipeline already has robust stochastic data flow:
stage-wise `update_frac`, `sampled` and `updatecnt` bookkeeping,
probabilistic pre-alignment, exact subset reproduction, online partial
reconstruction, and trailing reconstruction. That is not the missing piece.

The missing piece is an explicit optimizer contract:

- define the stochastic objective being optimized
- expose the mini-batch reconstruction as a stochastic normal-equation update
- separate particle sampling fraction from optimizer learning rate
- add optimizer state, learning-rate schedules, optional momentum, and
  diagnostics
- keep pose updates initially as stochastic latent-variable assignment through
  the existing polar matcher

The first implementation target is therefore a native SIMPLE
preconditioned-SGD volume optimizer, not an external autodiff engine and not a
replacement for the matcher.

## Current Machinery To Preserve

### Existing stochastic update flow

The current `abinitio3D` workflow already provides the stochastic outer loop:

1. `simple_commanders_abinitio.f90` derives the stage/run particle target from
   `nsample`, `nstates`, selected particles, and `UPDATE_FRAC_MAX`.
2. `simple_abinitio_controller.f90` emits per-stage `refine3D` command lines,
   including `update_frac`, `fillin`, `frac_best`, `balance`, `trail_rec`, and
   the active search mode.
3. `simple_matcher_smpl_and_lplims.f90` chooses the active particle subset.
4. `simple_oris_sampling.f90` records the subset through `sampled` and
   `updatecnt`.
5. `simple_commanders_prob.f90` samples once in `prob_align` and downstream
   `prob_tab`/`refine3D_exec` reproduce the exact same subset.
6. `simple_strategy3D_matcher.f90` updates orientations/states/shifts and
   writes partition-local current partial reconstructions.
7. `simple_commanders_rec_distr.f90` reduces partials, restores even/odd
   volumes, computes FSCs, postprocesses, and applies trailing reconstruction.

This flow should remain intact.

### Existing numerical kernels

The first SGD implementation must continue to use:

- `simple_reconstructor.f90`
  - `insert_plane_oversamp_opt`
  - `rho` / `rho_exp`
  - `sampl_dens_correct`
  - Kaiser-Bessel windowing through `simple_kbinterpol`
- `simple_reconstructor_eo.f90`
  - even/odd reconstruction composition
  - FSC and sampling-density correction
  - previous half-map/rho compatibility
- `simple_polarft_calc` and related `pftc` modules
  - polar reference storage
  - memoized objective evaluation
  - fast in-plane and shift search
  - existing analytic shift gradients
- `simple_matcher_refvol_utils.f90`
  - materialization of even/odd reprojection model files
  - active polar range handoff through reprojection model headers

## Mathematical Contract

For a mini-batch `B` and current assigned or sampled latent variables
`theta_i`, `shift_i`, and `state_i`, define the batch objective:

```text
L_B(V) = 1/2 sum_{i in B} || A_i V - y_i ||^2 / sigma_i^2 + R(V)
```

where:

- `V` is a state half-map in Fourier/volume representation
- `A_i` is the current projection, CTF, shift, symmetry, and interpolation
  operator for particle `i`
- `y_i` is the particle image
- `sigma_i^2` is the noise model used by the Euclidean objective
- `R(V)` represents optional regularization terms already expressed through
  ML regularization, FSC-derived priors, or nonuniform filtering policy

The current KB reconstruction path already computes the batch normal-equation
terms:

```text
b_B   = A_B^* y_B
rho_B ~= diag(A_B^* A_B)
V_B   = rho_B^{-1} b_B
```

The preconditioned stochastic gradient is:

```text
g_B(V_t) = rho_B V_t - b_B
```

A diagonal-preconditioned SGD step is:

```text
V_{t+1} = V_t - eta_t rho_B^{-1} g_B(V_t)
        = V_t + eta_t (V_B - V_t)
```

This is the key implementation bridge. The current partial reconstruction plus
sampling-density correction already produces `V_B`. Current trailing
reconstruction resembles the final line, but uses the realized update fraction
as a blending factor. To make the machinery bona fide SGD, the blend weight
must become an optimizer learning rate `eta_t`, distinct from `update_frac`.

## Core Design Decision

Keep two quantities separate:

```text
update_frac  = which particles are sampled
sgd_eta      = how far the map moves toward the mini-batch target
```

The nominal or realized update fraction is sampling policy. It must remain
owned by the existing ab initio/refine machinery.

The learning rate is optimizer policy. It should be owned by a small
volume-domain SGD optimizer object or by explicit strategy/assembly state.

## User-Facing Parameter Surface

Every new command-line key must have a typed field in `simple_parameters`.

Proposed initial parameters:

```text
sgd               yes/no       default no
sgd_eta           real         default 0.2
sgd_eta_min       real         default 0.02
sgd_eta_decay     const/cos/inv_sqrt/stage
sgd_momentum      real         default 0.0
sgd_precond       batch/running/shell
sgd_pose_mode     assign/sample/soft/local
sgd_diag          yes/no       default yes
```

Recommended MVP subset:

```text
sgd               yes/no
sgd_eta           real
sgd_eta_decay     const/cos/inv_sqrt
sgd_diag          yes/no
```

Do not overload:

- `update_frac`
- `ufrac_trec`
- `frac_best`
- `trail_rec`

These already have workflow meanings and should remain compatible with current
ab initio and refine behavior.

## Ownership Boundaries

### Stage policy

Owner:

- `src/main/simple_abinitio_controller.f90`
- `src/main/commanders/simple/simple_commanders_abinitio.f90`

Responsibilities:

- decide whether a stage enables SGD
- pass `sgd=yes` and schedule parameters to child `refine3D`
- keep stage search policy unchanged unless deliberately testing an SGD variant
- preserve existing `update_frac`, `fillin`, `frac_best`, and `trail_rec`
  decisions

### Particle-domain update

Owner:

- `src/main/strategies/search/simple_strategy3D_matcher.f90`
- `src/main/strategies/search/simple_strategy3D_*.f90`
- `src/main/commanders/simple/simple_commanders_prob.f90`
- `src/main/simple_eul_prob_tab*.f90`

Responsibilities:

- choose or consume the sampled subset
- run polar matching/probability-table assignment
- update particle orientations/states/shifts
- preserve single-read batch reconstruction
- write current mini-batch partial reconstruction artifacts

The particle-domain path should not know about `sgd_eta` in the MVP.

### Volume-domain SGD update

Owner:

- `src/main/commanders/simple/simple_commanders_rec_distr.f90`
- new helper module, likely under `src/main/volume/`

Responsibilities:

- reduce current partials into mini-batch even/odd reconstruction targets
- restore `V_B_even` and `V_B_odd`
- read previous `V_t_even` and `V_t_odd` when SGD is enabled
- apply `V_{t+1} = V_t + eta_t (V_B - V_t)`
- write updated even/odd half-maps and merged state map
- emit SGD diagnostics

### Numerical reconstruction kernels

Owner:

- `src/main/volume/simple_reconstructor.f90`
- `src/main/volume/simple_reconstructor_eo.f90`
- `src/main/interp/simple_kbinterpol.f90`
- `src/main/interp/simple_gridding.f90`

Responsibilities:

- preserve the current KB splat/gather implementation
- expose additional low-level state only through narrow APIs if raw-gradient
  mode is later implemented
- avoid duplicating interpolation logic in commanders or strategies

## Proposed New Module

Add a small optimizer module:

```text
src/main/volume/simple_vol_sgd_optimizer.f90
```

Candidate public type:

```text
type :: vol_sgd_optimizer
contains
    procedure :: new
    procedure :: eta
    procedure :: apply_halfmap_update
    procedure :: write_diag
    procedure :: kill
end type
```

MVP behavior:

```text
V_new = V_prev + eta * (V_batch - V_prev)
```

This module should not read project files, choose particle subsets, or know
about probability tables. It should operate on `image` objects or narrowly
defined volume-domain inputs supplied by `volassemble`.

Future extensions can add:

- momentum half-maps
- shell-wise update normalization
- running diagonal preconditioner
- Adam/RMSProp-style second moments
- stochastic objective proxy calculations

## Refactor Plan

### Phase 0: Documentation and compile-time guardrails

Add this implementation note and then add a short policy update only after the
first code path exists.

Checklist:

- no behavior changes
- identify current calls where `trail_rec` does map blending
- identify all call paths that require previous half-maps

### Phase 1: Parameter plumbing

Add typed fields and parser registration:

- `src/main/params/simple_parameters.f90`
- `src/main/params/simple_parameters_parse.f90`
- `src/main/params/simple_parameters_phases.f90`
- UI registration if desired in `src/main/ui/simple/simple_ui_refine3D.f90`

Initial fields:

```text
character(len=3) :: sgd = 'no'
character(len=16) :: sgd_eta_decay = 'const'
real :: sgd_eta = 0.2
real :: sgd_eta_min = 0.02
logical :: l_sgd = .false.
logical :: l_sgd_diag = .true.
```

Validation:

- `sgd_eta > 0`
- `sgd_eta <= 1` for MVP convex-update mode
- `sgd_eta_min >= 0`
- `sgd_eta_min <= sgd_eta`
- only allow `sgd=yes` with `volrec=yes`
- initially require `trail_rec=yes` or explicit previous half-map availability

### Phase 2: Optimizer helper

Create `simple_vol_sgd_optimizer.f90`.

MVP API:

```text
call sgd_opt%new(params, state, which_iter)
eta = sgd_opt%eta()
call sgd_opt%apply_halfmap_update(prev_even, batch_even, new_even)
call sgd_opt%apply_halfmap_update(prev_odd,  batch_odd,  new_odd)
call sgd_opt%write_diag(...)
```

The first version can be stateless except for parameters and iteration index.
That keeps the initial integration simple and makes `sgd_eta=realized_fraction`
easy to compare against legacy trailing reconstruction.

### Phase 3: Split map blending from trailing reconstruction

Current relevant routine:

- `simple_commanders_rec_distr.f90::restore_state_from_parts`

Current conceptual sequence:

```text
reduce partials
restore even/odd current batch targets
restore merged current batch target
read restored halfmaps into build%vol/build%vol2
trail_restored_halves_if_needed
postprocess
```

Refactor toward:

```text
reduce partials
restore even/odd current batch targets
restore merged current batch target
read restored halfmaps into build%vol/build%vol2
apply_volume_feedback_if_needed
postprocess
```

where `apply_volume_feedback_if_needed` dispatches:

```text
if params%l_sgd:
    apply_sgd_volume_update
else:
    trail_restored_halves_if_needed
```

Important: this is a behavior-changing refactor. Keep the first code change
small and preserve the old path when `sgd=no`.

### Phase 4: Implement MVP SGD update in volassemble

When `sgd=yes`:

1. Require previous even/odd half-maps.
2. Read previous half-maps using the same compatibility logic already used by
   trailing reconstruction.
3. Treat the restored current half-maps as `V_B_even` and `V_B_odd`.
4. Compute:

```text
V_new_even = V_prev_even + eta_t * (V_B_even - V_prev_even)
V_new_odd  = V_prev_odd  + eta_t * (V_B_odd  - V_prev_odd)
```

5. Write:

```text
refine3D_state_halfvol_fname(state, 'even')
refine3D_state_halfvol_fname(state, 'odd')
refine3D_state_vol_fname(state)
```

6. Preserve nonuniform auxiliary source behavior by applying the same update to
   the NU source half-map images that currently receive trailing updates.

Comparison mode:

```text
sgd=yes sgd_eta=<realized_update_frac>
```

should closely reproduce the legacy trailing blend.

### Phase 5: Diagnostics

Add diagnostics to the existing benchmark/report style.

Per state, per iteration:

```text
SGD enabled
sgd_eta
update_frac command target
realized update fraction
||V_B_even - V_prev_even||
||V_B_odd  - V_prev_odd||
||eta * (V_B_even - V_prev_even)||
||eta * (V_B_odd  - V_prev_odd)||
relative update norm
state population
FSC 0.5 / 0.143 after update
```

Output options:

- append to existing `refine3D_volassemble_bench_fname`
- add `sgd_stateNN_iterMM.txt` only if detailed diagnostics are enabled

MVP should avoid storing large diagnostic volumes by default.

### Phase 6: Learning-rate schedules

Implement `sgd_eta_decay`.

Initial schedules:

```text
const:
    eta_t = sgd_eta

inv_sqrt:
    eta_t = max(sgd_eta_min, sgd_eta / sqrt(t_stage))

cos:
    eta_t = sgd_eta_min + 0.5 * (sgd_eta - sgd_eta_min) *
            (1 + cos(pi * progress))
```

Use existing decay helper patterns where practical, but keep optimizer policy
localized.

The stage controller may reset optimizer time at stage boundaries. That should
be explicit:

- `t_stage` for per-stage schedule
- `t_global` for full command invocation schedule

MVP recommendation: use per-stage `t_stage`.

### Phase 7: Optional momentum

Add only after MVP diagnostics are stable.

Momentum update:

```text
delta_t = V_B - V_t
m_t     = beta * m_{t-1} + (1 - beta) * delta_t
V_{t+1} = V_t + eta_t * m_t
```

State storage options:

1. in-memory only for shared-memory run
2. on-disk per state/even-odd momentum volumes for distributed/resumable runs

Do not add momentum to the first MVP unless the on-disk state contract is clear.

### Phase 8: Soft/top-K stochastic gradient

The MVP uses hard or sampled latent assignments. That is already a valid
stochastic generalized-EM/preconditioned-SGD hybrid.

For a stronger Bayesian-style stochastic gradient, extend probability table
output to preserve top-K compact responsibilities:

```text
particle i:
    candidate k:
        state
        projection
        in-plane
        shift
        weight
```

Then reconstruct:

```text
b_B   = sum_i sum_k w_ik A_ik^* y_i
rho_B = sum_i sum_k w_ik diag(A_ik^* A_ik)
```

This is the natural extension of existing `prob` / `prob_neigh` candidate
machinery and avoids full posterior storage.

Implementation implications:

- `simple_eul_prob_tab*.f90` must emit compact top-K candidate assignments
- `simple_strategy3D_matcher.f90` or a new reconstruction helper must be able
  to backproject one particle into multiple candidate poses with weights
- the online single-read contract must still hold

This phase is intentionally after the hard-assignment SGD MVP.

### Phase 9: Continuous pose refinement

Do not start here.

Initial pose updates should continue to use the fast polar matcher and
probabilistic assignment machinery.

Later options:

1. Use existing analytic shift gradients in `simple_pftc_shsrch_grad`.
2. Fit local orientation tangent updates from polar-score neighborhoods around
   the best projection direction.
3. Add analytic derivatives through the KB projection operator for true
   continuous viewing-direction gradients.

The third option is a deeper research implementation. It should not be coupled
to the first map-SGD milestone.

## Detailed Data Flow For MVP

### Before matcher

Existing flow:

```text
materialize_reprojection_model
prob_align if enabled
read_reprojection_model
memoize_refs
```

No SGD changes.

### During matcher

Existing flow:

```text
build_batch_particles3D
choose_and_run_strategy
update project orientation/state/shift
prep_imgs4rec
update_rec
write_partial_recs
```

No SGD changes in MVP.

The current partials are interpreted as the mini-batch sufficient statistics
for the volume update.

### During volassemble

Existing flow:

```text
reduce partition-local partials
restore current even/odd halfmaps
restore current merged map
optionally trail with previous halfmaps
postprocess
```

SGD flow:

```text
reduce partition-local partials
restore current even/odd halfmap targets V_B
read previous halfmaps V_t
apply V_{t+1} = V_t + eta_t (V_B - V_t)
write updated even/odd halfmaps
write updated merged map
postprocess
```

### Next iteration

Existing flow:

```text
materialize_reprojection_model from current Cartesian halfmaps
write even/odd polar reprojection model files
matcher consumes polar references
```

No SGD changes.

## Compatibility Requirements

### `sgd=no`

Must reproduce current behavior.

### `sgd=yes, sgd_eta=1`

Should write the restored mini-batch target directly:

```text
V_{t+1} = V_B
```

This is intentionally aggressive and useful mainly for testing.

### `sgd=yes, sgd_eta=realized_update_frac`

Should match current trailing reconstruction within numerical tolerance,
modulo any ordering differences introduced by the refactor.

### `sgd=yes, update_frac=1`

Becomes full-batch preconditioned gradient/fixed-point update:

```text
V_{t+1} = V_t + eta_t (V_full - V_t)
```

This is useful for convergence diagnostics and high-resolution polishing.

## Tests And Validation

### Unit-level tests

Add targeted tests for `vol_sgd_optimizer`:

- `eta` schedule values
- `eta` bounds
- convex image update:

```text
prev = 10
batch = 20
eta = 0.25
new = 12.5
```

- `eta=0` rejected or produces no-op only if explicitly allowed
- `eta=1` returns the batch target

### Integration tests

Use a small synthetic or existing toy project:

1. Run current `abinitio3D` or `refine3D` path with `sgd=no`.
2. Run `sgd=yes, sgd_eta=realized_update_frac` and compare halfmaps after
   `volassemble`.
3. Run `sgd=yes, sgd_eta=1` and confirm output equals mini-batch restored
   target.
4. Verify:
   - no second particle-stack read
   - orientation table updated exactly once per active particle
   - `sampled` and `updatecnt` semantics unchanged
   - reprojection model headers remain compatible
   - even/odd half-map names and FSC files are still produced

### Scientific diagnostics

Track:

- FSC after each update
- state-local update fraction
- map update norm per shell
- pose change distribution
- objective proxy, if implemented
- convergence rate against current abinitio3D baseline

The first success criterion is not immediate resolution gain. It is that the
explicit optimizer reproduces current behavior under compatible settings and
then gives stable, controllable behavior as `sgd_eta` varies.

## Risks

### Confusing sampling and optimization

Risk:

- treating `update_frac` as a learning rate

Mitigation:

- keep `update_frac` as subset selection
- introduce `sgd_eta` as the only map learning-rate field
- diagnostics must print both

### Breaking trailing reconstruction behavior

Risk:

- refactoring `restore_state_from_parts` changes legacy `trail_rec`

Mitigation:

- preserve the old path exactly under `sgd=no`
- add comparison tests using `sgd_eta=realized_update_frac`

### Overfitting mini-batch noise

Risk:

- large `eta` chases noisy mini-batch maps

Mitigation:

- conservative default `sgd_eta`
- decay schedules
- optional shell-wise update diagnostics
- preserve even/odd separation

### Introducing hidden optimizer state

Risk:

- momentum or running preconditioner becomes ambient global state

Mitigation:

- keep optimizer state builder-owned or explicit
- serialize state only when resumability is deliberately implemented

### Disrupting polar matching performance

Risk:

- inserting SGD logic into matcher or reference preparation slows the hot path

Mitigation:

- MVP does not touch polar matching kernels
- reference handoff remains unchanged
- SGD is volume-domain assembly logic

## Recommended Milestone Order

1. Add parameter plumbing and `sgd=no` no-op behavior.
2. Add `vol_sgd_optimizer` with stateless convex half-map update.
3. Refactor `volassemble` to dispatch legacy trailing versus SGD update.
4. Add diagnostics and comparison mode.
5. Validate `sgd_eta=realized_update_frac` against current trailing behavior.
6. Tune default `sgd_eta` and schedules on small datasets.
7. Add optional momentum only after the MVP is stable.
8. Add compact top-K soft stochastic reconstruction.
9. Investigate continuous orientation-gradient refinement.

## Initial MVP Definition

The first mergeable version should support:

```text
sgd=yes
sgd_eta=<0..1>
sgd_eta_decay=const
sgd_diag=yes/no
```

and implement:

```text
V_new_even = V_prev_even + eta * (V_batch_even - V_prev_even)
V_new_odd  = V_prev_odd  + eta * (V_batch_odd  - V_prev_odd)
```

It should not implement:

- momentum
- Adam/RMSProp
- top-K soft reconstruction
- continuous orientation gradients
- external autodiff

That restraint is deliberate. It proves the optimizer contract while preserving
the mature abinitio3D stochastic sampling, KB reconstruction, and polar matching
machinery.

