# Continuous 3-D PCG post-reconstruction pose-polishing PLAN

**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md)

## How to enable and use the feature

The feature is disabled by default. Run it after refinement by enabling one
joint rotation-and-shift pass in `reconstruct3D`:

```bash
simple_exec \
  prg=reconstruct3D \
  projfile=/path/to/project.simple \
  rec_backend=pcg \
  pcg_pose_polish=yes \
  combine_eo=no \
  pgrp=d2 \
  mskdiam=180 \
  nthr=8
```

The option is valid only for `reconstruct3D` with an active PCG backend and
independent even/odd maps. SIMPLE first reconstructs fixed PCG half-maps from
the completed project's poses. It then refines each active particle against
only its matching half-map, persists only accepted complete poses, and runs one
additional two-iteration PCG reconstruction. No orientation-search iteration
runs. Internal workers do not start another polish.

The log must contain one `PCG POSE POLISH` terminal summary and one
`RUNNING FINAL PCG RECONSTRUCTION` marker. Omit the option, or use
`pcg_pose_polish=no`, to keep the established behavior.

For the complete Oracle production validation, use
[`production/tests/continuous_3D_pcg_pose_validation/README.md`](../../production/tests/continuous_3D_pcg_pose_validation/README.md).
The production-contract runner stores one default audit arm, four matched
beta-gal A/B arms, and all evidence in one timestamped directory. The separate
truth-diagnostic runner described below stores its complete clean/noisy and
exact/perturbed matrix in one timestamped directory.

## Current decision: production gate failed

Oracle run `continuous_3D_pose_validation_20260815_120459` passed all ten
numerical test groups, all four invalid-command checks, default-off behavior,
and shared/distributed equivalence. It did not pass the scientific A/B gate:

- all 2,000 particles were accepted as improved;
- the reported particle objective decreased from `19.8838` to `19.3898`;
- the shared and distributed FSC-area change was `-0.0112817`;
- cFAR decreased from `0.8495` to `0.8399`;
- FSC=0.143 resolution changed from `3.671` to `3.744` Angstrom;
- the initial exported poses matched the simulation truth within the
  orientation-table export precision.

This is a real, reproducible acceptance failure, not a validator formatting
error. The optimizer reduced its same-particle objective while it moved exact
poses away from truth and reduced independent half-map agreement. The result
does not by itself identify a wrong derivative or disprove local pose
polishing. It leaves two primary hypotheses: mismatch between the independent
simulation operator and the PCG polishing operator, and high-frequency noise
overfitting. Both can be present.

Do not weaken the FSC tolerance. Keep `pcg_pose_polish` disabled by default and
do not recommend production use until the truth-diagnostic matrix below passes.

## Implementation goal

Own `pcg_pose_polish=yes` in `reconstruct3D`. Reuse the validated joint
five-parameter mathematics, particle I/O, half-map ownership, weighting, and
persistence. Remove the separate `refine3D` terminal hook so reconstruction is
the single production owner.

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

Status terms are evidence-specific: `SOURCE COMPLETE` means lightweight source
checks passed, `ORACLE PASS` means the supplied Oracle output was observed, and
`SCIENTIFIC FAIL` means the executable completed but a predeclared scientific
criterion failed. Line numbers identify the current source checkpoint and must
be refreshed after later edits.

Repository checkpoints group several phases because commits were created only
after their Oracle gates:

| Checkpoint | Scope |
| --- | --- |
| `222160351` | Test workbox, fixtures, fast-KB derivative, packed gather, and PCG workspace foundations |
| `686090c73` | Expanded synthetic half-set, conventional-control, trajectory, and lambda validation |
| `52de96f5c` | Shift/joint pose numerical and initial production-polisher implementation |
| `1e8be6002` | Shift, rotation, recovery, and batch tests |
| `2e00ea5b9` | FINAL joint-pose SPEC and initial PLAN |
| current worktree | `reconstruct3D` ownership correction, pose contract, production A/B package, and truth-diagnostic matrix |

### Phase 1 — isolated mother/child test workbox

**Goal:** create one mother executable with independently compiled child cases,
while keeping shared helpers free of substantive test logic.

**Major file changes:**

- Mother runner, `case=` selection, child execution, and suite accounting:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:1` and `:60`.
- Shared assertions, deterministic seed, and reusable truth builder:
  `production/tests/simple_continuous_3D_pcg_refinement_helpers.f90:27`, `:59`,
  and `:76`.
- First executable contract:
  `production/tests/simple_continuous_3D_pcg_refinement_scaffold_test.f90:10`.
- Dedicated child-module ownership was established for volume, noise, half-set,
  KB, shift, rotation, and recovery cases. Their numerical bodies were filled
  in later phases. CMake already globbed these sources, so no CMake file changed.

**Progress and result:** Oracle Work Package 1 passed on 2026-08-11. The mother
ran eight groups: one pass, seven explicit status-77 skips, and zero failures.
Log SHA-256: `C7687C6F8E86BD272A0A61BE90A19DBD84438D2097C3120E257F027F00619771`.

**Status:** `ORACLE PASS`. This proved the test architecture, not any numerical
method.

### Phase 2 — deterministic volume and noise fixtures

**Goal:** produce a repeatable asymmetric 3-D phantom, add controlled noise,
and verify independent noisy observations before testing reconstruction.

**Major file changes:**

- Truth-volume builder shared by later cases:
  `production/tests/simple_continuous_3D_pcg_refinement_helpers.f90:76`.
- Volume dimensions, fingerprints, variance, norm, and asymmetry:
  `production/tests/simple_continuous_3D_pcg_refinement_volume_test.f90:15`.
- Noise case coordinator:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_test.f90:11`.
- Unit-Gaussian replacement:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_gauran_test.f90:18`.
- Added volume noise and replay checks:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_volume_test.f90:21`.
- Independent noisy production-projector observations:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_observation_test.f90:33`.
- Shared statistical calculations only:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_support.f90:1`.

**Progress and result:** Work Package 2 passed with the expected deterministic
volume fingerprints. The first Work Package 3 build failed because three
`image%set_rmat` calls omitted the required `ft` argument; therefore the old
binary's skips were not accepted as evidence. After passing `.false.` for the
real-space inputs, Oracle Work Package 3 passed. Requested SNR 0.5 realized as
`0.501052/0.494972`; requested observation SNR 1.0 realized as
`1.009640/0.982848`; independent-noise correlations remained close to zero.
Final log SHA-256: `8A583349A3E951D57559612C5E7FE11D363E6FD4F70E89E25CEB4042AA1354CC`.

**Status:** `ORACLE PASS`. The independent observations were a fixture gate,
not yet a PCG forward-model oracle.

### Phase 3 — derivative preflight, fast KB, and packed gather

**Goal:** differentiate the functions that SIMPLE actually executes and define
a fixed-volume workspace that takes one padded Fourier snapshot.

**Major file changes:**

- Executed degree-14 polynomial value and derivative:
  `src/main/interp/simple_kbinterpol.f90:181`.
- Normalized 3-D fixed-cell stencil and its three gradients:
  `src/main/interp/simple_kbinterpol.f90:343`.
- One-snapshot Fourier workspace:
  `src/main/volume/simple_reconstructor_pcg.f90:931`.
- Packed/Friedel value and gradient gather:
  `src/main/volume/simple_reconstructor_pcg.f90:989`.
- Polynomial, normalized stencil, support, and switch tests:
  `production/tests/simple_continuous_3D_pcg_refinement_kb_test.f90:20`.
- Packed gather and differentiated Friedel tests:
  `production/tests/simple_continuous_3D_pcg_refinement_kb_gather_test.f90:19`.

**Progress and result:** the design preflight fixed the `rho` terminology,
right-tangent and shift-phase conventions, fixed-cell derivative scope,
gauge/symmetry handling, weak-particle outcomes, and LM policy before numerical
implementation. Fast KB/stencil tests passed with log SHA-256
`F36D8524A0766CC981D13E1086BCF1594A5BF0C40EA88B5F064B7698E2ABAC42`.
The packed gather then passed with value relative error `0`, derivative finite-
difference error `2.332225e-7`, and differentiated Friedel error
`1.870813e-9`. Packed-gather log SHA-256:
`616853ECABB64B97EA8A3EF07C7B39F3D21823ACBF5DE79BB3A38576AD077CF7`.

**Status:** `ORACLE PASS`. These tests establish local derivatives inside one
fixed stencil cell; they do not establish complete pose recovery.

### Phase 4 — Stage 0 independent half-set reconstruction

**Goal:** validate half ownership and characterize PCG semi-convergence and
regularization with observations generated by the production projector rather
than by the PCG gather under test.

**Major file changes:**

- Half-set test coordinator and reported trajectories:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_test.f90:35`.
- Disjoint orientations, independent observations, PCG trajectories, fixed
  support, residuals, FSC, and MRC output:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_support.f90:35`,
  `:63`, `:144`, `:169`, `:200`, and `:272`.
- Conventional gridding control on identical observations:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_gridding.f90:21`.
- Consolidated view-count, iteration, lambda, residual, norm, truth-error, FSC,
  and volume-output matrix:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90:36`.

**Progress and result:** the first source had invalid Fortran continuations in
seven assertion calls; compilation passed after those test-only corrections.
Early noisy truth-correlation failures led to noiseless controls, fixed-support
controls, conventional gridding, forced iteration trajectories, and a bracketed
lambda sweep. The data showed PCG semi-convergence: useful noisy recovery
peaked near 4–8 iterations and degraded by iteration 40. The final scale-
sensitive sweep selected lambda `3000`; noisy raw errors were
`0.562039/0.555794`, compared with conventional `0.690137/0.679725`.
Final log SHA-256: `DA40358B35F9D533DC7EEC787B32672C96E606CCE36C1BC5E470C49D0DE0B3E1`.

**Status:** synthetic Stage 0 is `ORACLE PASS`. It established half isolation
and exposed stopping/regularization sensitivity. It did not prove that
standalone PCG is superior, and it did not validate pose polishing.

### Phase 5 — Stage 1 shift derivative and bounded LM

**Goal:** validate the two shift coordinates first, then exercise them through
a half-isolated batch scaffold and production weighting.

**Major file changes:**

- Shift residual, Jv, real adjoint, normal terms, and objective gradient:
  `src/main/volume/simple_reconstructor_pcg.f90:1012`, `:1044`, `:1072`,
  `:1106`, and `:1148`.
- Bounded two-parameter LM and explicit outcomes:
  `src/main/volume/simple_reconstructor_pcg.f90:1297`.
- Production transfer and observation whitening used by the tests:
  `src/main/volume/simple_reconstructor_pcg.f90:1717` and `:1786`.
- Fixed-half shift batch and atomic rollback:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:50`.
- Shift sign, adjoint, gradient, recovery, exact, weak, and weighting tests:
  `production/tests/simple_continuous_3D_pcg_refinement_shift_test.f90:23`.
- Batch half-isolation and terminal-accounting test:
  `production/tests/simple_continuous_3D_pcg_refinement_shift_polish_test.f90:17`.

**Progress and result:** the CTF-free/unit-noise scaffold passed first. The
first production-weighted build then failed at a `tiny(1.0)` sigma guard because
SIMPLE's `TINY` symbol hid the intrinsic; the finite-positive sigma check was
corrected. Oracle then passed `shift_gradient`, `shift_polish`, and the mother
suite. Weighted gradient/recovery/whitening errors were
`1.632803e-10/1.848185e-6/2.396121e-8`; six particles improved, two exact
solutions remained unchanged, and the largest trial step was `0.2765353` pixel.
Log SHA-256: `C782301FD4F71E19FC87EF7882100F3B288D8AB3DA264FB8962812016D061076`.

**Status:** `ORACLE PASS` as component evidence. No public shift-only feature
remains because rotation and shift errors are coupled.

### Phase 6 — joint rotation-and-shift numerics

**Goal:** replace the shift-only production concept with one five-parameter
right-tangent pose update and one scaled bounded LM solve.

**Major file changes:**

- Right-increment rotation composition:
  `src/main/volume/simple_reconstructor_pcg.f90:289`.
- Rotation Jv, joint normal terms, objective gradient, and stencil telemetry:
  `src/main/volume/simple_reconstructor_pcg.f90:1167`, `:1202`, `:1256`, and
  `:1276`.
- Scaled bounded five-parameter LM:
  `src/main/volume/simple_reconstructor_pcg.f90:1430`.
- Joint batch update and complete-pose rollback:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:127`.
- Rotation finite differences and independent right-increment oracle:
  `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90:21`.
- Joint recovery, exact-solution, weak-system, and rollback cases:
  `production/tests/simple_continuous_3D_pcg_refinement_recovery_test.f90:20`.
- Mother-suite registration:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:60`.
- Frozen contract and this implementation plan:
  `doc/implementation_notes/continuous_3D_pose_end_polishing_spec.md:1` and
  `doc/implementation_notes/continuous_3D_pose_end_polishing_plan.md:1`.

**Progress and result:** the first Oracle compile found three `tiny(...)` calls;
they were replaced by a double-precision numerical floor. The numerical gate
then passed. Joint recovery reduced rotation/shift error from
`0.03148287` rad/`0.4920366` pixel to `0`/`3.202989e-7` pixel. Review found one
production persistence defect: rejected poses were still rewritten through an
Euler-matrix round trip. Persistence now occurs only for
`POSE_LM_ACCEPTED_IMPROVEMENT`, and the rebuilt suite passed. Final log SHA-256:
`1C311AFFA4CD32A252B91ED0C43D5391177ACE2D8C384771DAF86A574CC57E03`.

**Status:** joint numerical core is `ORACLE PASS`.

### Phase 7 — production contract and ownership

**Goal:** make `reconstruct3D` the only production owner, keep the option off by
default, protect worker routes, persist only accepted poses, and run one final
PCG reconstruction.

**Major file changes:**

- Default parameter storage:
  `src/main/params/simple_parameters.f90:308`.
- Command-line registration:
  `src/main/params/simple_parameters_parse.f90:152`.
- Typed activation and rejection rules:
  `src/main/params/simple_parameters_phases.f90:801`.
- `reconstruct3D` UI exposure and PCG-only activation:
  `src/main/ui/simple/simple_ui_refine3D.f90:138`.
- Base reconstruction, final polish, project persistence, and second PCG pass:
  `src/main/commanders/simple/simple_commanders_rec.f90:33`.
- Removal of the earlier terminal `refine3D` ownership attempt:
  `src/main/commanders/simple/simple_commanders_refine3D.f90:1336`.
- Distributed worker stripping so only the master polishes:
  `src/main/strategies/parallelization/simple_rec3D_strategy.f90:307`.
- Removal of obsolete `refine3D` child/master propagation:
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90:143` and
  `:1077`.
- Production fixed-half particle loading, CTF/whitening, accepted-pose write,
  normalization-mask fallback, terminal summary, and active Fourier limit:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:214`, `:268`, `:306`,
  `:337`, and `:355`.
- Point-group, global-gauge, multi-orientation, half-isolation, bounds, and
  terminal-accounting contract test:
  `production/tests/simple_continuous_3D_pcg_refinement_pose_contract_test.f90:22`.
- Ten-case mother-suite registration:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:60` and `:183`.

**Progress and result:** the first production runner used terminal `refine3D`
arms. That design produced unmatched stochastic controls and exposed an
unallocated master normalization mask. Ownership was moved to deterministic
`reconstruct3D`; the master now creates a fallback normalization mask when
`build%lmsk` is absent. Oracle run `continuous_3D_pose_validation_20260815_120459`
passed all ten numerical groups, all invalid-command checks, default-versus-
explicit-off behavior, persistence/consumption checks, and shared/distributed
equivalence.

**Status:** production software contract is `ORACLE PASS`; scientific utility
is not.

### Phase 8 — frozen beta-gal production A/B

**Goal:** compare matched `pcg_pose_polish=no|yes` arms from one frozen noisy
beta-gal checkpoint without running an orientation-search iteration.

**Major file changes:**

- One-directory, collect-all Oracle runner, policy cases, evidence export, and
  shared/distributed arms:
  `production/tests/continuous_3D_pcg_pose_validation/run_oracle_validation.sh:157`,
  `:172`, `:216`, and `:228`.
- Pose/metadata/FSC parsing, A/B checks, default behavior, and route equivalence:
  `production/tests/continuous_3D_pcg_pose_validation/analyze_pose_ab.py:96`,
  `:211`, `:277`, and `:298`.
- Reproducible command, fixed tolerances, artifact layout, and interpretation:
  `production/tests/continuous_3D_pcg_pose_validation/README.md:19`, `:48`,
  `:69`, and `:93`.

**Progress and result:** the first policy runner omitted required arguments and
mistook a UI-input stop for a typed-policy result. It was corrected to verify
that each invalid case reaches the intended validation layer and to continue
all independent work after a failure. The corrected run completed. All 2,000
particles were accepted and the internal objective decreased from `19.8838` to
`19.3898`, but FSC area changed by `-0.0112817` on both routes, cFAR changed
from `0.8495` to `0.8399`, and FSC=0.143 resolution changed from `3.671` to
`3.744` Angstrom. The starting poses matched the simulation truth within export
precision.

**Status:** `SCIENTIFIC FAIL`. This is not a harness failure. The current
same-particle acceptance objective can move truth poses and reduce independent
map agreement. Keep the feature disabled by default.

### Phase 9 — clean/noisy truth-diagnostic matrix

**Goal:** distinguish an operator/objective consistency defect, local pose-
solver failure, and noisy high-frequency overfitting before changing production
numerics.

**Major file changes:**

- Deterministic perturbation of only `e1`, `e3`, `x`, and `y` while preserving
  CTF and half metadata:
  `production/tests/continuous_3D_pcg_pose_validation/prepare_truth_oritab.py:32`.
- Deterministic alternating half ownership when the simulator table has no
  split, plus microscope-metadata validation:
  `production/tests/continuous_3D_pcg_pose_validation/prepare_truth_oritab.py:46`
  and `:57`.
- One collect-all clean/noisy, exact/perturbed, full/FSC-limited 16-arm runner;
  source manifest; clean production-projector stack; project creation; FSC
  evidence; and final packaging:
  `production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh:132`,
  `:158`, `:234`, `:272`, and `:293`.
- Truth-pose errors, half-map and truth-map FSC changes, per-case decisions, and
  aggregate feature gate:
  `production/tests/continuous_3D_pcg_pose_validation/analyze_truth_matrix.py:34`,
  `:52`, `:95`, and `:186`.
- User procedure and interpretation:
  `production/tests/continuous_3D_pcg_pose_validation/README.md:109`.

The eight matched cases are:

| Observations | Initial poses | Polish band | Question |
| --- | --- | --- | --- |
| clean | exact truth | full | Is the complete pipeline stationary at independently generated truth? |
| clean | perturbed | full | Can the local optimizer recover known pose errors without noise? |
| noisy | exact truth | full | Does full-band polishing overfit noise? |
| noisy | perturbed | full | Does full-band polishing recover useful pose signal? |
| noisy | exact truth | FSC=0.5 | Does a conservative band preserve truth poses? |
| noisy | perturbed | FSC=0.5 | Does that band retain useful recovery? |
| noisy | exact truth | FSC=0.143 | Does a wider supported band preserve truth poses? |
| noisy | perturbed | FSC=0.143 | Does that wider band recover useful signal? |

Predeclared limits are `0.001745` rad exact-pose rotation RMS increase,
`0.05` pixel exact-pose shift RMS increase, and `0.005` maximum FSC-area
decline. Clean perturbed rotation and shift errors must both decrease by at
least 50 percent. Noisy perturbed errors must both decrease by at least 10
percent. The aggregate gate requires both clean cases and at least one complete
noisy-band pair to pass. Do not tune these limits after reading the result.

**Progress and result:** the first Oracle run
`continuous_3D_pose_truth_diagnostic_20260817_092117` did not execute a
scientific arm. All four `import_particles` commands stopped during UI input
validation because `kv`, `cs`, and `fraca` were not supplied, but SIMPLE
returned status zero. Empty projects then reached `partition_eo` and failed at
index zero; the missing `poses_after.txt` files were downstream consequences.
The runner now treats UI-input text as failure, supplies and validates the
microscope values, assigns deterministic alternating half ownership when the
input table has only `eo=0`, and requires the imported particle count and both
half counts to pass before scheduling any arm.

**Status:** corrected `SOURCE COMPLETE`; Bash syntax, Python AST, perturbation
metadata, whitespace, and conflict-marker checks passed. The first Oracle run
is a harness failure and provides no scientific evidence. Oracle rerun is
pending. This is the current phase.

## Current next command

The complete command, inputs, output directory, and result-reading order are in
[`production/tests/continuous_3D_pcg_pose_validation/README.md`](../../production/tests/continuous_3D_pcg_pose_validation/README.md#truth-diagnostic-matrix-after-an-ab-failure).
Run only `run_truth_diagnostic.sh`; it schedules every Phase 9 arm and keeps all
evidence in one timestamped directory.

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

The component implementation and production connection have passed their
focused, mother-suite, persistence, and shared/distributed gates. The feature
is not scientifically complete because the frozen beta-gal A/B failed. The
truth-diagnostic matrix must identify and resolve that failure before the SPEC
can be accepted. Stage 4 alternating reconstruction remains blocked.
