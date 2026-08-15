# Continuous 3-D joint pose end-polishing PLAN

**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md)

## Implementation goal

Replace the uncommitted shift-only public prototype with
`pcg_pose_polish=yes`. Reuse its validated shift mathematics, particle I/O,
half-map ownership, weighting, persistence, and final-reconstruction path while
adding the three rotation coordinates and one joint five-parameter solve.

## Validated prerequisites

- Shift residual, sign, Jacobian, adjoint, bounded LM, CTF weighting, whitening,
  and recovery passed on Oracle Linux.
- Fixed-half batch isolation and rollback passed on Oracle Linux.
- The executed fast KB derivative, normalized stencil gradient, packed gather,
  and stencil-switch measurements pass the existing `kb_derivative` case.
- Latest numerical log SHA-256:
  `C782301FD4F71E19FC87EF7882100F3B288D8AB3DA264FB8962812016D061076`.

These are component gates. The real-data production gate is intentionally
deferred until the complete joint pose system exists.

## Mathematical implementation

With row-vector sampling

$$
\ell_0=\mathrm{padf}\,[h,k,0]R,
$$

use a right tangent increment. At zero increment,

$$
\frac{\partial\ell}{\partial\omega_a}=\ell_0\times e_a,
$$

and the three rotation columns are

$$
J_{\omega_a}=T_iS_i\left(\nabla\widehat V(\ell_0)\cdot
(\ell_0\times e_a)\right).
$$

The two validated shift columns remain

$$
J_{t_x}=i\frac{2\pi h}{N}m_i,\qquad
J_{t_y}=i\frac{2\pi k}{N}m_i.
$$

Accumulate the real five-by-five Gauss--Newton block and gradient in double
precision. Solve in dimensionless coordinates. One rotation scale is the angle
that moves the mask radius by one cropped pixel,

$$
s_\omega=\frac{1}{\max(r_{\mathrm{mask}},1)}\ \mathrm{rad},\qquad
s_t=1\ \mathrm{pixel}.
$$

Use these scales for damping and step bounds. Retain the validated gain-ratio
policy: accept only positive actual reduction with ratio at least `0.25`, reduce
damping above `0.75`, and otherwise increase damping. Every trial recomputes the
full objective and stencil.

## Work packages

Current implementation and validation status:

- Work packages 1--3 and the shared production connection in work package 5
  are implemented.
- `rotation_gradient` and the C1 joint-recovery, exact, weak, bound, rollback,
  and telemetry parts of work package 4 pass on Oracle Linux.
- The focused rotation run and the complete nine-case mother suite passed on
  Oracle Linux with no skipped or failed groups. Log SHA-256:
  `1C311AFFA4CD32A252B91ED0C43D5391177ACE2D8C384771DAF86A574CC57E03`.
- Point-group recovery, explicit global-gauge fixtures, half-isolation for
  rotations, distributed equivalence, persistence, and beta-gal remain open.
- A source-review correction now prevents rejected or unreliable poses from
  being rewritten through an Euler-matrix round trip. The complete Oracle
  Linux suite passed after this correction.

### 1. Public-contract rename

- Replace `pcg_shift_polish` with `pcg_pose_polish` in parameters, parsing, UI,
  validation, child-command stripping, commander activation, and diagnostics.
- Rename the production strategy and its summary to joint-pose terminology.
- Keep shift-only numerical routines and tests as internal component coverage.

### 2. Rotation derivative API

- Extend `pcg_fourier_workspace` with weighted rotation and joint-pose normal
  terms using its executed packed gather gradient.
- Add a stable right-increment $SO(3)$ update and re-linearize after every
  accepted step.
- Preserve the active Fourier limit and record stencil-switch proximity/counts.

### 3. Joint LM

- Assemble one scaled symmetric five-by-five block per particle.
- Add finite checks, eigenvalue or rank diagnostics, diagonal floors, angular
  and shift bounds, gain-ratio acceptance, and explicit terminal outcomes.
- Roll back both rotation and shift unless the terminal result contains an
  accepted improvement.

### 4. Focused tests

- Implement `case=rotation_gradient` with fixed-cell centered differences and
  one-sided checks near stencil switches.
- Implement `case=pose_recovery` with independent known rotations and shifts,
  exact/no-information cases, gauge fixing, and point-group-aware errors.
- Extend batch tests to prove rotation/shift immutability across halves and
  balanced outcome accounting.

### 5. Production integration

- Generalize the current final polisher to read and write complete poses.
- Persist Euler metadata through existing orientation/project APIs without
  changing state, half labels, or discrete-search decisions.
- Run the existing additional two-iteration PCG reconstruction with accepted
  poses, once on the master only.

### 6. Oracle and production gates

Run:

```bash
simple_test_continuous_3D_pcg_refinement case=rotation_gradient \
  2>&1 | tee continuous_3D_pcg_refinement.log

simple_test_continuous_3D_pcg_refinement case=pose_recovery \
  2>&1 | tee -a continuous_3D_pcg_refinement.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_refinement.log
```

After numerical acceptance, clone one frozen beta-gal checkpoint. Change only
`pcg_pose_polish=no|yes`, first on the shared route and then on the distributed
route. Compare persisted poses, terminal counts, final-reconstruction evidence,
angular/shift distributions, and independent FSC.

## Review and verification roles

- The implementation pass owns source changes and lightweight non-compiling
  checks.
- A fresh review pass checks derivative signs, $SO(3)$ convention, scaling,
  symmetry/gauge treatment, half ownership, and production orchestration.
- The user runs Oracle compilation, focused tests, the mother suite, and frozen
  beta-gal workflows.

## Risks

- A wrong left/right increment can pass self-generated recovery tests; finite
  differences must perturb the executed rotation matrix independently.
- Centered differences that cross a KB stencil switch do not test a local
  derivative.
- Poor radians/pixels scaling can make a correct five-by-five system unusable.
- Symmetry-equivalent poses and global gauge can create false recovery failures.
- Rotation/shift coupling can allow a lower objective with a physically wrong
  pose; report both coordinates and image-space displacement.

## Completion gate

The joint pose implementation is complete only after focused and mother tests
pass on Oracle Linux, production poses persist into the final PCG reconstruction,
shared/distributed behavior agrees, and the frozen beta-gal A/B satisfies the
SPEC. Stage 4 alternating reconstruction remains blocked until then.
