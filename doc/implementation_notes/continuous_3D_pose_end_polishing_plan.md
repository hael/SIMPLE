# Continuous 3-D PCG post-reconstruction pose-polishing PLAN

**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md)

## Implementation standing and Hans decision boundary

The isolated implementation is software-complete for its current FINAL SPEC:
the option is default-off, the production path compiles, the numerical and
operator contracts passed, serial/parallel pose results agree, and a clean
64-thread production matrix completed. The new 500-particle amortized batch
also compiles; its Oracle runtime timing remains pending.

The method is not scientifically accepted. Phase 10F failed exact-pose
rotation stationarity, clean joint recovery, and noisy FSC gates. This blocks
production recommendation and normal-workflow integration, but it does not
block a focused commit of the disabled implementation and its evidence.

Hans does not need to approve the current code snapshot. Input is required only
to authorize a next design:

1. **Workflow:** stop here, retain a one-off terminal command, or redesign it
   as a selected `abinitio3D`/`refine3D_auto` stage?
2. **Reference:** continue with the reconstructed same-half fixed map, or design
   a statistically independent/regularized reference contract?
3. **Frequency and rotation:** keep the current native-Nyquist PCG operator, or
   specify one shared reconstruction-and-polishing cutoff or rotation prior?
4. **Acceptance:** is a lower particle objective sufficient, or must updates
   pass a held-out, FSC-neutral, or other dataset-level criterion?

No new numerical policy should be implemented until those answers produce a
new or revised FINAL SPEC.

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
  nthr=64
```

The option is valid only for `reconstruct3D` with an active PCG backend and
independent even/odd maps. SIMPLE first reconstructs fixed PCG half-maps from
the completed project's poses. It then refines each active particle against
only its matching half-map, persists only accepted complete poses, and runs one
additional two-iteration PCG reconstruction. No orientation-search iteration
runs. Internal workers do not start another polish.

This is an isolated top-level command contract, not current SIMPLE workflow
integration. The standard `SIMPLE_data_testing` scripts call `abinitio3D` and
`refine3D_auto`; they do not call `prg=reconstruct3D` directly or propagate
`pcg_pose_polish`. Therefore the present implementation does not polish after
the internal reconstructions of those workflows. Whether polishing belongs
once at the end or after selected reconstruction iterations is a product and
scientific decision for Hans and the other refinement owners.

The log must contain one `PCG POSE POLISH` terminal summary and one
`RUNNING FINAL PCG RECONSTRUCTION` marker. Omit the option, or use
`pcg_pose_polish=no`, to keep the established behavior.

For the active Oracle scientific gate, use the truth-matrix command in the
[pose-validation README](../../production/tests/continuous_3D_pcg_pose_validation/README.md).
It stores the complete clean/noisy and exact/perturbed matrix plus the frozen
contract, plan, source hashes, fixture hashes, and handoff in one timestamped
directory. The separate beta-gal production-contract runner is retained but
deferred because it exercises a direct command rather than current workflow
integration.

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

### Operator-contract review after Phase 10B

The source and document review on 2026-08-17 found two distinct mismatches that
must not be combined under one `reference bias` label:

1. **Confirmed PCG-policy mismatch.** The PCG normal operator predicts from
   the inverse-envelope volume, `G(E^-1 V)`, but the production pose polisher
   passes the reconstructed half-map directly to `set_volume`. This violates
   the original plan's explicit deapodization requirement and the current PCG
   operator policy. The mismatch is visible at
   `src/main/volume/simple_reconstructor_pcg.f90:1906` and
   `src/main/strategies/search/simple_pcg_pose_polisher.f90:347`.
2. **Confirmed validation-generator mismatch.** `simulate_particles` creates
   a padded projection and clips it in real space before the native FFT. The
   PCG gather omits that finite-box operation by design. The existing
   `matched_window` control proves `W G(V)` for the clean simulator fixture; it
   does not prove that `W` belongs in the production PCG likelihood.
3. **Unresolved production-model choice.** Conventional SIMPLE reconstruction
   and PFTC matching use finite observed particle boxes but do not explicitly
   pass each trial reference through the simulator's hard-clip operator. The
   simulator-matched `W` is therefore not yet established as the standard
   experimental-data likelihood.
4. **Incomplete production preprocessing match.** The clean matched-window
   control does not include the complete observed-particle normalization and
   data-dependent edge statistics. A production correction must freeze any
   statistics measured from the observed particle and reuse them for the
   complete batch. The executed reconstruction input path is
   `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90:299`; the
   corresponding polishing input path starts at
   `src/main/strategies/search/simple_pcg_pose_polisher.f90:380`.
5. **Frequency-policy mismatch remains possible.** The PCG reconstruction
   operator is initialized at native Nyquist, while the pose workspace applies
   `params%kfromto` afterward. Both must use one authoritative shell setting.

The FINAL SPEC now separates the fixed data preprocessing
$\bar y_i=Q_i(I_i;\eta_i)$ from the production PCG prediction operator. With
$V_h=P u_h$ denoting the already supported stored half-map,

$$
\mathcal A_i(R_i,t_i)V_h=T_i(t_i)G_P(R_i)E^{-1}V_h,
\qquad T_i=C_iS_i/\sqrt{\sigma_i^2}.
$$

It requires operation order, inverse envelope, support, transfer, whitening,
and shell range to match in reconstruction, prediction, and every Jacobian
column. Data-derived preprocessing statistics stay fixed. The simulator's
finite-box operator remains a validation hypothesis: adding it only to the
polisher is forbidden because that would violate the same-operator contract.
Selecting it requires a coordinated forward, adjoint, normal-operator, and
pose-Jacobian redesign.

#### Reconciliation of the four design documents

- The frozen original proposal and `doc/policies/reconstruct3D_pcg_policy.md`
  agree on the current production contract: the solved volume is supported and
  inverse-envelope corrected before the PCG gather; CTF, shift, and whitening
  then act in the packed Fourier plane.
- The FINAL SPEC previously said only "same operator" and did not state whether
  particle-box preprocessing belonged on the data side or model side. The
  clarification above closes that ambiguity and forbids a polishing-only
  finite-box transform.
- `continuous_3D_pose_end_polishing_Literature Review and Recommended Path.md`
  gives a principled finite-box operator $W$ as an alternative likelihood. It
  is a design candidate, not evidence that current SIMPLE reconstruction
  already executes $W$.
- The present source follows the PCG data-input path for normalization, taper,
  FFT, extraction, CTF, and whitening, but omits $E^{-1}$ when it snapshots the
  half-map for polishing. That omission is the immediate P1 implementation
  mismatch. Whether the PCG likelihood itself needs $W$ is a separate P2 design
  question that Phase 10C measures before discussion with Hans.

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
| `1c7d5191b` | Move the opt-in production owner from `refine3D` to `reconstruct3D` |
| `50f26f9c6` | Pose contract, production A/B package, and truth-diagnostic matrix |
| `2c167c557` | Nine-phase implementation and evidence record |
| `3ad305e19` | Production-routine and important numerical-operation documentation |
| `fd1f5d24e` | Corrected Phase 9 scientific result and next isolation experiment |
| current worktree | Phase 10D operator correction, final truth-matrix package, workflow-boundary finding, and scientific handoff |

### Phase 1 — isolated mother/child test workbox

**Goal:** create one mother executable with independently compiled child cases,
while keeping shared helpers free of substantive test logic.

**Major file changes:**

- Mother runner, `case=` selection, child execution, and suite accounting:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:1` and `:63`.
- Shared assertions, deterministic seed, and reusable truth builder:
  `production/tests/simple_continuous_3D_pcg_refinement_helpers.f90:28`, `:63`,
  and `:81`.
- First executable contract:
  `production/tests/simple_continuous_3D_pcg_refinement_scaffold_test.f90:11`.
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
  `production/tests/simple_continuous_3D_pcg_refinement_helpers.f90:81`.
- Volume dimensions, fingerprints, variance, norm, and asymmetry:
  `production/tests/simple_continuous_3D_pcg_refinement_volume_test.f90:16`.
- Noise case coordinator:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_test.f90:12`.
- Unit-Gaussian replacement:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_gauran_test.f90:19`.
- Added volume noise and replay checks:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_volume_test.f90:22`.
- Independent noisy production-projector observations:
  `production/tests/simple_continuous_3D_pcg_refinement_noise_observation_test.f90:34`.
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
  `src/main/volume/simple_reconstructor_pcg.f90:992`.
- Polynomial, normalized stencil, support, and switch tests:
  `production/tests/simple_continuous_3D_pcg_refinement_kb_test.f90:21`.
- Packed gather and differentiated Friedel tests:
  `production/tests/simple_continuous_3D_pcg_refinement_kb_gather_test.f90:20`.

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
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_test.f90:36`.
- Disjoint orientations, independent observations, PCG trajectories, fixed
  support, residuals, FSC, and MRC output:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_support.f90:36`,
  `:65`, `:148`, `:174`, `:206`, and `:283`.
- Conventional gridding control on identical observations:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_gridding.f90:22`.
- Consolidated view-count, iteration, lambda, residual, norm, truth-error, FSC,
  and volume-output matrix:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90:37`.

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
  `src/main/volume/simple_reconstructor_pcg.f90:1015`, `:1047`, `:1075`,
  `:1109`, and `:1153`.
- Bounded two-parameter LM and explicit outcomes:
  `src/main/volume/simple_reconstructor_pcg.f90:1304`.
- Production transfer and observation whitening used by the tests:
  `src/main/volume/simple_reconstructor_pcg.f90:1726` and `:1795`.
- Fixed-half shift batch and atomic rollback:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:53`.
- Shift sign, adjoint, gradient, recovery, exact, weak, and weighting tests:
  `production/tests/simple_continuous_3D_pcg_refinement_shift_test.f90:24`.
- Batch half-isolation and terminal-accounting test:
  `production/tests/simple_continuous_3D_pcg_refinement_shift_polish_test.f90:18`.

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
  `src/main/volume/simple_reconstructor_pcg.f90:1172`, `:1207`, `:1263`, and
  `:1283`.
- Scaled bounded five-parameter LM:
  `src/main/volume/simple_reconstructor_pcg.f90:1437`.
- Joint batch update and complete-pose rollback:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:132`.
- Rotation finite differences and independent right-increment oracle:
  `production/tests/simple_continuous_3D_pcg_refinement_rotation_test.f90:22`.
- Joint recovery, exact-solution, weak-system, and rollback cases:
  `production/tests/simple_continuous_3D_pcg_refinement_recovery_test.f90:21`.
- Mother-suite registration:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:63`.
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
  `src/main/commanders/simple/simple_commanders_rec.f90:36`.
- Removal of the earlier terminal `refine3D` ownership attempt:
  `src/main/commanders/simple/simple_commanders_refine3D.f90:1336`.
- Distributed worker stripping so only the master polishes:
  `src/main/strategies/parallelization/simple_rec3D_strategy.f90:307`.
- Removal of obsolete `refine3D` child/master propagation:
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90:143` and
  `:1077`.
- Production fixed-half particle loading, CTF/whitening, accepted-pose write,
  normalization-mask fallback, terminal summary, and active Fourier limit:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:244`, `:282`, `:294`,
  `:302`, `:380`, `:401`, and `:427`.
- Point-group, global-gauge, multi-orientation, half-isolation, bounds, and
  terminal-accounting contract test:
  `production/tests/simple_continuous_3D_pcg_refinement_pose_contract_test.f90:23`.
- Ten-case mother-suite registration:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:63` and `:193`.

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
  `production/tests/continuous_3D_pcg_pose_validation/run_oracle_validation.sh:159`,
  `:175`, `:220`, and `:233`.
- Pose/metadata/FSC parsing, A/B checks, default behavior, and route equivalence:
  `production/tests/continuous_3D_pcg_pose_validation/analyze_pose_ab.py:104`,
  `:227`, `:294`, and `:316`.
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
  `production/tests/continuous_3D_pcg_pose_validation/prepare_truth_oritab.py:34`.
- Deterministic alternating half ownership when the simulator table has no
  split, plus microscope-metadata validation:
  `production/tests/continuous_3D_pcg_pose_validation/prepare_truth_oritab.py:49`
  and `:61`.
- One collect-all clean/noisy, exact/perturbed, full/FSC-limited 16-arm runner;
  source manifest; clean production-projector stack; project creation; FSC
  evidence; and final packaging:
  `production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh:135`,
  `:163`, `:240`, `:279`, and `:301`.
- Truth-pose errors, half-map and truth-map FSC changes, per-case decisions, and
  aggregate feature gate:
  `production/tests/continuous_3D_pcg_pose_validation/analyze_truth_matrix.py:34`,
  `:53`, `:99`, and `:192`.
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

The corrected Oracle run `continuous_3D_pose_truth_diagnostic_20260817_093049`
completed all 16 reconstruction arms and reached the analyzer. The analyzer was
the only recorded failure. The result is therefore scientific evidence, not a
harness failure.

The aggregate gate failed:

- `clean_exact_full` moved all 2,000 truth-start particles. Rotation RMS changed
  from `8.66466e-07` to `0.0171622` rad, shift RMS changed from `4.06222e-05`
  to `0.105219` pixels, half-map FSC area changed by `-0.00980282`, and
  truth-average FSC area changed by `-0.0211268`.
- `clean_perturbed_full` reduced shift RMS from `0.320157` to `0.106914` pixels
  but increased rotation RMS from `0.0173295` to `0.0223992` rad. It therefore
  failed the required joint 50-percent recovery gate.
- Full-band noisy exact and perturbed cases failed. The exact case moved truth
  poses and reduced half-map FSC area by `-0.0112817`.
- The FSC=0.5 and FSC=0.143 perturbed cases passed and improved FSC, but both
  matched exact-pose controls failed their stationarity tolerances.
- Every enabled arm reported 2,000 improved particles and zero unchanged,
  unreliable, step-bound, invalid, or iteration-limit outcomes. A lower local
  objective therefore did not identify whether a pose change was physically
  correct.

**Status:** `SCIENTIFIC FAIL`. The clean exact result supports a mismatch
between the independent truth generator and the fixed-map polishing objective.
Possible causes include an operator/objective inconsistency and bias in the
reconstructed half-map used as the fixed reference. The current matrix does not
separate those causes. Full-band noisy results also support noise overfitting.
Keep the feature disabled by default.

#### Why the corrected Phase 9 run failed

The 16 reconstruction arms completed normally. Only the scientific analyzer
failed. The enabled command first made two-iteration PCG half-maps, froze each
matching half-map, minimized the CTF-weighted fixed-map Fourier residual for
each particle, persisted accepted poses, and reconstructed again. LM accepted
a trial when the fully recomputed fixed-map objective decreased and its gain
ratio was at least `0.25`. The acceptance rule did not use the known truth pose
or FSC.

The `clean_exact_full` arm is decisive. The runner used `snr=10`; SIMPLE's
simulator adds Gaussian noise only when `snr < 5`, so this was a no-noise
control. All 2,000 particles started at truth within export precision, all
2,000 were classified as improved, and the aggregate objective decreased from
`277.299` to `268.286`. The physical result moved in the opposite direction:
rotation RMS reached `0.0171622` rad, shift RMS reached `0.105219` pixels,
half-map FSC area decreased by `0.00980282`, and truth-average FSC area
decreased by `0.0211268`. Therefore, lower residual against the reconstructed
fixed half-map is not a valid proxy for a correct pose in this configuration.

The current evidence supports these hypotheses but does not yet choose between
them:

1. **Reconstruction-map bias.** The reference is an approximate, regularized,
   cropped two-iteration PCG half-map, not the volume that generated the clean
   observations. Pose changes can compensate for errors in that map.
2. **Leave-in bias.** Each particle contributes to the same half-map against
   which it is polished. Gold-standard half isolation prevents cross-half
   leakage but does not remove particle-to-own-map correlation.
3. **Forward-path mismatch.** The clean images use SIMPLE's production
   projector. Polishing uses the PCG packed Fourier gather, fast KB
   interpolation, crop, CTF, whitening, phase, and precision conventions.
   Self-consistent derivative tests do not prove equality to the independent
   generator.
4. **Rotation-specific bias or coupling.** In `clean_perturbed_full`, shift RMS
   improved from `0.320157` to `0.106914` pixels while rotation RMS worsened
   from `0.0173295` to `0.0223992` rad.
5. **Full-band noise overfitting.** Frequency-limited perturbed noisy cases
   improved pose error and FSC, while the full-band noisy exact case drifted
   and lost FSC. The frequency-limited exact controls still failed, so a cutoff
   alone is not yet an accepted correction.

Stencil switches remain a numerical risk, but they do not explain the result
by themselves. Some frequency-limited recovery arms passed despite millions
of switches, while the clean full-band arm failed with far fewer switches.

#### Frozen SPEC audit after Phase 9

The FINAL SPEC remains correct and is unchanged. Acceptance criterion 2
requires exact poses to remain unchanged. Acceptance criterion 8 requires
neutral or improved independent half-map FSC. Phase 9 violated both criteria;
it did not invalidate them. The failure does not require a new public option or
a changed production contract. Any later proposal to change the objective,
reference ownership, or acceptance criterion requires discussion with Hans and
formal SPEC change control before production implementation resumes.

### Phase 10 — fixed-reference operator isolation

**Goal:** locate the first source of truth-pose drift by separating the KB
envelope, independent forward generator, LM/Jacobian path, and reconstructed
reference.

**Major file changes:**

- Eleven-arm collect-all fixed-reference matrix and diagnosis ladder:
  `production/tests/simple_continuous_3D_pcg_refinement_fixed_reference_test.f90:47`.
- Simulator-matched padded/projected/clipped generator and independent padded
  exact-KB generator:
  `production/tests/simple_continuous_3D_pcg_refinement_fixed_reference_test.f90:227`
  and `:258`.
- Shared Stage 0 half-set observations now use the same padded projection and
  real-space clipping convention as `simulate_particles`:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_support.f90:65`.
- PCG-matched control, shared measurement/polishing path, amplitude-fitted
  residuals, and low/mid/high-shell evidence:
  `production/tests/simple_continuous_3D_pcg_refinement_fixed_reference_test.f90:295`,
  `:312`, `:337`, and `:411`.
- Explicit diagnostic-only child selector, excluded from the mother suite
  while the failure is unresolved:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:197`.
- Three-stage simulator boundary diagnostic, transition metrics, and diagnosis
  ladder:
  `production/tests/simple_continuous_3D_pcg_refinement_forward_path_test.f90:31`,
  `:102`, `:150`, `:165`, and `:247`.
- Reused exact-pose measurement API and explicit `case=forward_path` selector:
  `production/tests/simple_continuous_3D_pcg_refinement_fixed_reference_test.f90:17`
  and `production/tests/simple_test_continuous_3D_pcg_refinement.f90:199`.
- Explicit matched-window model, transformed five-column Jacobian, derivative
  check, and test-only bounded recovery:
  `production/tests/simple_continuous_3D_pcg_refinement_matched_window_test.f90:61`,
  `:274`, `:430`, `:486`, and `:700`.
- Shared simulator-stage generator and diagnostic-only `case=matched_window`
  selector:
  `production/tests/simple_continuous_3D_pcg_refinement_forward_path_test.f90:16`
  and `production/tests/simple_test_continuous_3D_pcg_refinement.f90:202`.
- Commands, controlled variables, result interpretation, and limits:
  `production/tests/continuous_3D_pcg_pose_validation/README.md:182`, `:232`,
  and `:269`.

The diagnostic now compares three forward generators: the earlier native-box
helper, the padded-and-real-space-clipped `simulate_particles` geometry, and an
independent padded exact-KB Fourier slice. Raw truth,
inverse-envelope-corrected truth, globally amplitude-scaled raw and corrected
truth, PCG-matched corrected truth, and raw and corrected two-iteration
reconstructions form eleven arms. Each arm prints raw and amplitude-fitted residuals,
low/mid/high-shell residuals and correlations, pose-gradient norms, final pose
errors, objective change, and terminal outcomes.

Only numerical integrity, terminal accounting, and the already-established
PCG-matched LM control remain hard assertions. The open hypotheses are printed
in full before a deterministic diagnosis label. `EVIDENCE COMPLETE` means the
matrix ran; it is not a scientific acceptance result.

Run on Oracle Linux:

```bash
simple_test_continuous_3D_pcg_refinement case=fixed_reference \
  2>&1 | tee continuous_3D_pcg_fixed_reference.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_fixed_reference.log
```

The first command collects all eleven operator arms. The second command is the
single regression pass for the shared half-set generator plus the established
KB, gradient, recovery, pose-contract, and production-contract cases.

Interpret the result in this order:

- Matched corrected-truth drift identifies the LM/Jacobian path.
- Stable matched poses but drifting amplitude-scaled exact-KB poses identify
  the exact-versus-fast padded gather or packed-addressing boundary.
- Stable exact-KB poses followed by stable scaled-raw but drifting
  scaled-corrected simulator poses identify the inverse-envelope model choice.
- Drifting scaled-raw simulator poses identify padded projection, real-space
  clipping, or native-plane normalization.
- Stable scaled-raw truth poses but drifting raw reconstructed-reference poses
  enter reconstruction-map or leave-in bias.
- Low/mid/high-shell residuals locate the frequency band where the first
  operator disagreement appears.

The first Oracle version ran on 2026-08-17. Its raw truth arm drifted to
`0.1254565` rad rotation RMS and `0.3162917` pixel-coordinate RMS while reducing
the objective from `4.501951` to `4.479724`. The four-arm follow-up then showed:

- raw independent truth: residual `0.9394884`, amplitude scale `15.89514`, and
  final rotation/shift RMS `0.1254565` rad / `0.3162917` pixels;
- corrected independent truth: residual `0.9311228`, scale `13.52865`, and
  final RMS `0.1267380` rad / `0.3313934` pixels;
- PCG-matched corrected truth: zero residual and gradient, 48 unchanged poses,
  and `0.0001914519` rad / zero-pixel RMS;
- corrected reconstruction: residual `0.3578665`, scale `1.020554`, and final
  RMS `0.06590715` rad / `0.1455302` pixels.

This rejects the inverse-envelope-only explanation and confirms that the
Jacobian/LM is stationary for its executed PCG operator. It also exposed a test
methodology gap: the earlier helper projects a native `BOX^3` volume directly
to a native `BOX^2` plane, whereas `simulate_particles` pads the volume and
plane and then clips in real space. The expanded update tests that boundary and
global normalization in one run without modifying production.

The first expanded Oracle run did not reach any scientific arm.
`OSMPL_PAD_FAC*[BOX,BOX,1]` multiplied all three dimensions and constructed a
`48x48x2` padded particle image. `image%clip` rejected clipping it into the
native `24x24x1` image. The mother suite independently reproduced the same
fixture error in `halfset_fsc`; its other nine groups passed. The padded plane
is now explicitly `[OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1]`, matching the
existing conventional half-set path and keeping the image two-dimensional.

The corrected nine-arm Oracle run completed. The matched PCG arm was exactly
stationary. The scaled independent exact-KB arm was also stationary with fitted
residual `1.954106e-5`, which clears the LM/Jacobian and exact-versus-fast KB
boundaries. Simulator observations fit scaled raw truth better than scaled
corrected truth (`0.01378730` versus `0.06638209`), while the tested scaled
corrected arm drifted to `0.02174075` rad and `0.03665505` pixels. The corrected
reconstruction drifted to `0.07516358` rad and `0.2095203` pixels. Because the
matrix did not polish scaled raw truth or the raw reconstruction, its printed
`SIMULATOR_PADDING_OR_CLIP_MISMATCH` label is provisional. Those two missing
arms are now included.

The mother suite's only failure was the shared fixture's requested-SNR gate.
Padded noise was normalized before clipping, so its native-box realized SNR no
longer represented the declared `0.5` test condition. The clean forward path
remains padded and clipped, but white noise is again calibrated on the native
observations under test. This preserves the controlled Stage 0 SNR contract.

The clean-build eleven-arm Oracle run completed on 2026-08-17, and all ten
mother-suite groups passed. The scaled raw-truth simulator arm was far closer
to the PCG model than scaled corrected truth: fitted residual `0.01378730`
versus `0.06638209`. It nevertheless moved exact poses by `0.002032991` rad
and `0.003275161` pixels. The rotation error is slightly above the frozen
`0.001745` rad stationarity tolerance. Its fitted low/mid/high residuals were
`0.01001645`, `0.1222118`, and `0.9886179`, locating the disagreement mainly at
middle and high frequencies.

The raw and corrected two-iteration reconstruction arms showed larger bias.
They moved exact poses by `0.1044913` rad / `0.2146550` pixels and `0.07516358`
rad / `0.2095203` pixels, respectively. This confirms that inverse-envelope
correction is not the solution and that reconstructed-reference bias remains
after the earlier LM/Jacobian and exact-versus-fast KB controls passed.

The new `case=forward_path` diagnostic now removes LM and reconstruction. It
compares the raw-truth PCG gather with native frequencies captured directly
from the padded Fourier projection, after the simulator Fourier round trip,
and after real-space clipping/native FFT. Adjacent-stage fitted residuals will
identify the first executed boundary above the predeclared `1e-3` tolerance.

Run on Oracle Linux after rebuilding:

```bash
simple_test_continuous_3D_pcg_refinement case=forward_path \
  2>&1 | tee continuous_3D_pcg_forward_path.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_forward_path.log
```

The clean-build forward-path Oracle run completed on 2026-08-17. Direct padded
projection and Fourier-roundtrip transition residuals were `1.543602e-4` and
`1.281264e-4`, below the predeclared `1e-3` tolerance. The real-space clip and
native-FFT transition residual was `0.01378673`, with low/mid/high residuals
`0.01001685`, `0.1231536`, and `6.341939`. The deterministic diagnosis was
`REALSPACE_CLIP_NATIVE_FFT`. The appended mother suite passed all ten groups.

This establishes a discrete forward-model mismatch rather than an LM,
Jacobian, KB, or Fourier-roundtrip defect. Real-space clipping is a linear
window operation, so the test-only matched model applies the same
`F_native * clip * F_padded^-1` transformation to the padded PCG prediction and
all five pose derivatives. It then checks the transformed gradient, exact-pose
stationarity, and controlled joint recovery before any production design is
considered.

Run on Oracle Linux after rebuilding:

```bash
simple_test_continuous_3D_pcg_refinement case=matched_window \
  2>&1 | tee continuous_3D_pcg_matched_window.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_matched_window.log
```

The matched-window Oracle run completed on 2026-08-17. The fitted residual was
`1.345176e-6`; exact-pose rotation and shift errors were zero and
`1.369075e-7` pixels. The injected `0.01414014` rad / `0.25` pixel perturbation
recovered to zero rotation error and `2.451475e-7` pixels. The diagnosis was
`MATCHED_WINDOW_SUPPORTED`, and the appended mother suite passed all ten
groups. This confirms the clip/native-FFT correction for an exact truth
reference. It does not clear reconstructed-reference or leave-in bias.

**Status:** `MATCHED-WINDOW TRUTH CONTROL PASSED; REFERENCE-BIAS MATRIX PENDING`.
Production correction remains blocked.

### Phase 10B — reconstructed-reference and ownership isolation

**Goal:** determine why a reconstructed fixed half-map still moves exact poses
after the forward window is matched.

**Major file changes:**

- Reusable matched-window batch measurement with shell restriction, exact and
  perturbed recovery, and rigid-gauge removal:
  `production/tests/simple_continuous_3D_pcg_refinement_matched_window_test.f90:180`,
  `:274`, `:486`, and `:586`.
- Explicit single-half reconstruction helper with controlled lambda, iteration
  count, and support:
  `production/tests/simple_continuous_3D_pcg_refinement_halfset_support.f90:162`.
- Twenty-two-arm collect-all reference matrix and contribution diagnosis:
  `production/tests/simple_continuous_3D_pcg_refinement_reference_bias_test.f90:44`
  and `:185`.
- Diagnostic-only child selector:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:205`.
- One unattended Oracle runner with per-case logs and source/checksum manifests:
  `production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh:1`.
- Matrix contract, interpretation, and commands:
  `production/tests/continuous_3D_pcg_pose_validation/README.md:307`.

The matrix uses eight clean exact even-half probe particles. It compares the
truth volume, own-half and opposite-half references after 2/4/8 PCG
iterations, and a two-iteration reference reconstructed after excluding all
eight probes. It repeats the two-iteration comparison over low, middle, and
full shell ranges; masked and unmasked references; and lambda values `1e-3`,
`1`, `100`, and `2000`. All arms use the matched clipping operator.

Each arm reports amplitude-fitted residual, exact-pose gradients, raw pose
drift, drift after removal of the best common global rotation and projected
three-dimensional translation, perturbed-pose recovery, accepted steps, and
objective changes. Component flags identify a material ownership, gauge,
frequency, iteration, regularization, or support contribution. The independent
and holdout references have different angular coverage; agreement between
both is required before attributing an improvement specifically to leave-in
bias.

Only finite evidence and the already-passed truth-volume stationarity/recovery
controls are hard gates. Open scientific hypotheses remain collect-all.

Run all operator diagnostics and the mother suite in one directory:

```bash
bash production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh \
  --output-parent "$HOME/Projects"
```

Oracle package `continuous_3D_operator_diagnostics_20260817_143406` completed
with no harness failures. The fixed-reference, forward-path, matched-window,
reference-bias, and ten-case mother-suite markers were present. The
reference-bias matrix reported that iteration count and support materially
changed the drift, while gauge, ownership, frequency restriction, and lambda
did not meet its 25-percent component threshold. Two-iteration own-half drift
was `0.09357053` rad / `0.1902455` pixels; eight iterations reduced it to
`0.03736626` rad / `0.1372004` pixels. Removing support reduced the
two-iteration result to `0.06251357` rad / `0.09779254` pixels.

**Status:** `ORACLE EVIDENCE COMPLETE; INTERPRETATION SUPERSEDED BY PHASE 10C`.
The label `RECONSTRUCTED_REFERENCE_BIAS` is provisional because those arms
passed each reconstructed map directly to the gather and therefore still
included the confirmed missing inverse-envelope operation.

### Phase 10C — complete operator-contract matrix

**Goal:** separate the PCG inverse-envelope mismatch, simulator finite-box
mismatch, production edge taper, reconstruction iteration count, and synthetic
box/support margin before changing production numerics.

**Major file changes:**

- Dynamic independent projector, fixed-support truth, current matrix-free PCG
  reference reconstruction, four prediction models, and finite-difference
  stationarity measurements:
  `production/tests/simple_continuous_3D_pcg_refinement_operator_contract_support.f90:51`,
  `:60`, `:96`, `:150`, `:171`, and `:266`.
- Three-box, four-model, truth/reference collect-all matrix and decision flags:
  `production/tests/simple_continuous_3D_pcg_refinement_operator_contract_test.f90:27`
  and `:109`.
- Explicit `case=operator_contract` dispatch:
  `production/tests/simple_test_continuous_3D_pcg_refinement.f90:17` and `:208`.
- Unattended Oracle collection, checksum manifest, commands, interpretation,
  and limitations:
  `production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh:50`,
  `:72`, and `:79`; and
  `production/tests/continuous_3D_pcg_pose_validation/README.md:360`.

The declared models are

$$
G(V),\qquad G(E^{-1}V),\qquad WG(V),\qquad WG(E^{-1}V).
$$

The same supported asymmetric object is embedded in boxes 24, 32, and 48 while
its support radius remains 10 pixels. Each box uses 48 independent SIMPLE
projector observations. Exact-truth controls are measured before and after the
production edge taper. Tapered observations then produce 2-, 4-, and
8-iteration PCG references; each reference is evaluated under all four models.
Eight fixed probe orientations report fitted residual, nuisance amplitude,
normalized objective, and numerical rotation/shift gradient RMS.

The hard gates cover finite evidence and exact independent truth matching by
the explicit finite-box model. Reconstruction iteration and box-margin
comparisons remain collect-all evidence because they select the next production
design. Clean variance normalization is deliberately excluded:
noiseless solvent has no valid noise standard deviation. A later noisy
production diagnostic must measure normalization statistics from the observed
particle once and hold them fixed for every prediction and Jacobian evaluation.

Run the focused matrix or the complete unattended package:

```bash
simple_test_continuous_3D_pcg_refinement case=operator_contract \
  2>&1 | tee continuous_3D_pcg_operator_contract.log

bash production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh \
  --output-parent "$HOME/Projects"
```

Oracle package `continuous_3D_operator_diagnostics_20260817_180441` completed
with zero recorded failures. All five focused diagnostics emitted their
completion markers, all ten mother-suite cases passed, and no non-finite value,
warning, or floating-point exception was reported.

The finite-box truth control behaved as designed. Across boxes 24, 32, and 48,
$WG(V)$ matched the independently projected observations with fitted residuals
between `2.91e-7` and `6.10e-7`. The bare $G(V)$ residual decreased from
`5.06e-3` to `2.52e-3` as the object-to-edge margin increased. This confirms
the simulator clip boundary and shows that its effect shrinks with margin.

The deapodization flag from the reconstructed-reference matrix was false under
its predeclared 10-percent threshold. $WG(E^{-1}V)$ improved the reconstructed
reference residual by `5.68`--`8.51` percent in box 24, `2.29`--`5.29` percent
in box 32, and only `0.75`--`1.15` percent in box 48. This does not reject the
PCG inverse envelope. The independent generator calls
`projector%fproject_serial`, which evaluates the same normalized KB
interpolation family and therefore carries the same gather envelope. Its exact
$WG(V)$ match is an inverse-crime control for the envelope question. The PCG
policy explicitly says that only its envelope-free Stage 8 validates
deapodization.

The result supports three decisions:

1. Apply $E^{-1}$ exactly once to the already supported stored half-map before
   the polishing Fourier workspace. This is a confirmed same-PCG-operator P1,
   independent of the simulator matrix.
2. Preserve the preceding PCG solve's native-Nyquist shell range. The current
   pose-only `params%kfromto` override is a confirmed same-operator P2.
3. Do not add $W$ only to polishing. At realistic margin it has little effect
   on reconstructed-reference residual, and a polishing-only window would
   violate the FINAL SPEC. A future $W$ proposal still requires coordinated
   forward, adjoint, normal-operator, and pose-Jacobian work.

**Status:** `ORACLE MATRIX COMPLETE; P1/P2 SAME-OPERATOR CORRECTIONS SELECTED`.

### Phase 10D — PCG operator-contract corrections

**Goal:** make the production fixed-map pose objective use the same inverse
envelope and shell range as the PCG reconstruction, and add a non-projector
envelope control for the corrected objective.

**Major file changes:**

- Apply the inverse KB envelope once to the already supported half-map before
  `set_volume`; retain the workspace's reconstruction shell range instead of
  overriding it with `params%kfromto`:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:347` and `:350`.
- Build the established PCG Stage 8 envelope-free observation control
  $G(E^{-1}V)$:
  `production/tests/simple_continuous_3D_pcg_refinement_operator_contract_support.f90:151`.
- Require that the inverse-envelope pose model matches that control and that
  the bare gather does not:
  `production/tests/simple_continuous_3D_pcg_refinement_operator_contract_test.f90:29`
  and `:107`.
- Record why the independent SIMPLE projector is valid for the finite-box
  boundary but not independent for the KB envelope:
  `production/tests/continuous_3D_pcg_pose_validation/README.md:386`.

The envelope-free control deliberately reuses the already accepted
`test=pcg_recon` Stage 8 construction. It isolates the envelope contract; it
does not replace fixed-cell derivative tests or claim an independent KB
implementation. Production data normalization, CTF, whitening, half ownership,
LM policy, and the default-off activation contract are unchanged.

**2026-08-17 overnight validation plan:** synchronize this checkout to the
Oracle Linux validation checkout, rebuild the debug installation, and run the
complete unattended operator-diagnostic package. Keep all logs, manifests, and
the summary in its timestamped directory under `$HOME/Projects`.

This run has four acceptance gates:

1. the corrected production and diagnostic sources compile and link together;
2. the envelope-free control reports
   `clip/taper/envelope-control/margin T T T T`, proving that the executed
   inverse-envelope model matches the controlled observations while the bare
   gather does not;
3. the fixed-reference, forward-path, matched-window, reference-bias, and
   operator-contract cases complete so their evidence continues to separate
   KB-envelope, finite-box clip, taper, support-margin, iteration, and
   reconstructed-reference effects; and
4. the unchanged ten-case mother suite passes ten of ten.

Passing these gates establishes compilation, same-PCG-operator wiring, and
regression safety. It does not establish a scientific benefit from pose
polishing. After a pass, the next work is the truth-controlled production
matrix, the frozen simulated beta-gal matched A/B, and the raw-data HPC
beta-gal matched A/B, in that order. A failure stops that sequence and must be
diagnosed from the retained case log before another numerical change.

The 2026-08-17 UCRT64 preflight resolved the expected local source and rsync
executable. The first dry run found
`continuous_3D_matrix_volumes_20260817_180550614` and its MRC evidence inside
the remote validation checkout. That result directory was moved intact to its
correct location directly under `$HOME/Projects`. A second `--delete` dry run
had no deletion candidates, and the explicitly approved real synchronization
completed. Remote shell and Perl line endings were normalized. The first debug
install attempt compiled the changed Fortran but failed when an install-time
`simple_exec` check loaded the system C++ and Fortran runtimes. Repeating the
incremental build with
`LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64` matched the configured
GCC 15.2 compiler and completed installation.

Run after rebuilding on Oracle Linux:

```bash
bash production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh \
  --output-parent "$HOME/Projects"
```

The updated `operator_contract` case must report
`clip/taper/envelope-control/margin T T T T`; the mother suite must remain ten
of ten. Then rerun the truth-controlled production matrix and frozen beta-gal
A/B before enabling or recommending the feature.

Oracle package
`$HOME/Projects/continuous_3D_operator_diagnostics_20260817_192152` completed
in 325 seconds with zero recorded failures. Its canonical storage path is
`/usr/local/data/mazhar/Projects/continuous_3D_operator_diagnostics_20260817_192152`.
All result-manifest entries verified successfully. The focused cases emitted
their completion markers, the operator contract reported
`clip/taper/envelope-control/margin T T T T`, and the mother suite passed all
ten scheduled groups with zero skips or failures.

The complete package was also downloaded to
`C:\msys64\home\hossainm7\Projects\continuous_3D_operator_diagnostics_20260817_192152`.
Its local `MANIFEST.sha256` verification passed every retained input, log,
summary, and status file; `STATUS.txt` contains `PASS`.

The new envelope-free control separates the P1 wiring defect. For boxes 24,
32, and 48, the fitted residual for the inverse-envelope model was
`8.906049e-8`, `7.851385e-8`, and `7.433188e-8`; the bare gather residual was
`5.120360e-2`, `3.295628e-2`, and `1.235027e-2`. Thus the executed
$G(E^{-1}V)$ path matches the controlled observations and the old $G(V)$ path
does not.

The matrix also reported `reconstructed-reference deapod-improves-10pct F`.
This does not reverse the policy correction: applying $E^{-1}$ is required for
the same-PCG operator. It shows that deapodization alone does not remove the
larger reconstructed-reference mismatch. The earlier boundaries remain:
finite-box clip, taper, envelope control, and support margin are distinguishable;
the reference-bias matrix still flags iteration count and support, but not
gauge, leave-in ownership, frequency restriction, or lambda at its declared
25-percent threshold. No further production numerical change is selected from
this result.

**2026-08-18 plan:**

1. Run `run_truth_diagnostic.sh` with the frozen truth volume, truth
   orientations, noisy stack, and the previously measured 4.255/3.671-A FSC
   limits. Review `analysis/truth_matrix.md` and the clean/noisy,
   exact/perturbed, and full/FSC-limited arms. This is the first test of whether
   the corrected production operator keeps exact poses stationary and improves
   controlled perturbations.
2. Package the complete truth-controlled directory and use
   `continuous_3D_pose_end_polishing_handoff.md` for scientific review. Do not
   weaken a gate or make another numerical correction before that review.
3. Defer the frozen beta-gal and HPC explicit-command A/B runs. They can test
   the isolated `reconstruct3D` component, but they do not demonstrate normal
   `abinitio3D` or `refine3D_auto` workflow behavior because those workflows do
   not currently propagate `pcg_pose_polish`.

**Status:** `ORACLE OPERATOR CONTRACT AND TEN-CASE REGRESSION PASS`; the final
truth-controlled production matrix and scientific handoff remain pending.

### Phase 10E — particle-level pose-polish parallelism

**Goal:** make `nthr` control the expensive fixed-map particle LM work, not
only the two PCG reconstructions around it.

The first 64-thread truth-matrix run exposed a performance defect after the
base PCG solve. Process telemetry reported `%CPU=267`, which is approximately
2.67 occupied cores and can appear as about one percent of a large host. The
process retained approximately 15.7 GB because the project, particle images,
and reconstruction state remained resident. The stable memory footprint did
not indicate parallel progress. Source inspection confirmed that the pose LM
loop was serial and had no progress output until all 2,000 particles finished.

**Major file changes:**

- Parallelize independent per-particle LM systems against the immutable
  Fourier workspace, retain private solver state per particle, and reduce all
  summaries serially in particle order:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:158` and `:196`.
- Use the established `MAXIMGBATCHSZ=500` production image batch while the
  active OpenMP thread count controls simultaneous LM workers; report requested
  threads, active threads, and batch size:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:290`.
- Emit one progress record after each completed batch so a long arm can be
  distinguished from a stalled process:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:407`.
- Run identical mixed observable/unobservable batches with one thread and up
  to four threads, then require identical poses, terminal outcomes, counters,
  objective sums, and step bounds:
  `production/tests/simple_continuous_3D_pcg_refinement_recovery_test.f90:130`
  and `:179`.

For `nthr=64`, the expected production marker is:

```text
>>> PCG POSE POLISH THREADS: REQUESTED 64 ACTIVE 64 BATCH 500
```

The LM region should then use up to 64 cores when a half contains at least 64
particles. Each worker processes its static share of the 500-particle batch.
I/O, setup, summary reduction, pose persistence, and small final batches remain
partly or fully serial. The worker count is the active `nthr`; the batch size
does not create 500 OpenMP workers.

The running `continuous_3D_pose_truth_diagnostic_20260818_101659` process uses
both binaries because the Oracle executable was rebuilt while its unattended
runner was still scheduling arms. The early enabled arms were serial and took
`884.7` and `1108.5` seconds. Later enabled arms reported
`REQUESTED 64 ACTIVE 64 BATCH 64` and completed in `192.4`--`213.6` seconds.
The runner stopped during `noisy_perturbed_fsc05_on` after polishing and the
final two PCG half-solves. It contains only 12 of 16 arms, has no analysis, and
cannot support a scientific truth-matrix decision. It is retained only as
mixed-binary performance evidence.

The rebuilt Oracle `pose_recovery` case passed its one-thread versus
four-thread equivalence comparison. The complete mother suite then passed all
ten scheduled groups with zero skips and failures. This validates deterministic
parallel pose results and component regression safety. A clean 64-thread truth
matrix must now run without another rebuild or source synchronization. Compare
its scientific results with the serial arms and inspect the thread and progress
markers.

For a 200,000-particle production project, fixed-size batching keeps numerical
working memory independent of the total particle count. The source now uses
SIMPLE's established `MAXIMGBATCHSZ=500` batch while retaining `nthr`
concurrent LM workers. Progress follows SIMPLE's large-loop convention: one
record every five full batches and one at the end of each half. Oracle must
measure the new batch policy; do not infer large-N scaling from the
2,000-particle test.

**Status:** `64-PARTICLE-BATCH ORACLE COMPILATION, SERIAL/PARALLEL EQUIVALENCE,
AND TEN-CASE REGRESSION PASS; 500-PARTICLE BATCH SOURCE CHECKED; CLEAN TRUTH
MATRIX AND LARGE-N SCALING PENDING`.

### Phase 10F — clean parallel truth matrix

Oracle directory `continuous_3D_pose_truth_diagnostic_20260818_111338`
completed all 16 arms and produced the declared analysis. Every enabled arm
reported `REQUESTED 64 ACTIVE 64 BATCH 64` and completed in
`71.6`--`74.3` seconds; disabled controls completed in `19.4`--`20.2` seconds.
The one recorded failure is the scientific analyzer, not compilation,
integration, or harness execution.

The clean exact control improved substantially relative to the pre-Phase-10D
matrix but remained outside the frozen rotation tolerance. Exact rotation RMS
moved from `8.66e-7` to `0.0038695` rad and shift RMS moved from `4.06e-5` to
`0.0250316` pixels. The shift and FSC-area limits passed, but rotation exceeded
the `0.001745`-rad limit. The optimizer classified 1,998 particles as improved
and two at the iteration limit while reducing its fixed-map objective from
`576.585` to `573.949`. Thus the corrected fixed-map objective is still not
stationary at truth under the predeclared gate.

The clean perturbed control reduced rotation RMS only from `0.0173295` to
`0.0143828` rad, about 17 percent, while shift RMS improved from `0.320157` to
`0.0331244` pixels, about 90 percent. It failed the requirement that both
coordinates improve by at least 50 percent. Truth-average FSC area improved by
`0.0136901`; the failure is specifically rotation recovery.

The noisy exact full-band control drifted to `0.0035222` rad and `0.0468496`
pixels and reduced half-map/truth-average FSC area by `0.0130845`/`0.0223521`.
The noisy perturbed full-band control improved rotation by about 10 percent and
shift by about 45 percent, but truth-average FSC area declined by `0.00630986`,
slightly beyond the `0.005` tolerance.

The `full`, `fsc05`, and `fsc0143` noisy results are bit-for-bit identical in
pose errors, FSC changes, objectives, attempted and accepted steps, stencil
switches, and terminal outcomes. Their commands differ only by `lp`, showing
that those arms did not impose distinct shell limits on the executed PCG pose
objective. They cannot support a three-band frequency-policy conclusion. The
full-band controls remain valid evidence. A real frequency-policy experiment
would require one authoritative limit shared by reconstruction, prediction,
and every Jacobian column; a pose-only shell override remains forbidden by the
FINAL SPEC.

The feature gate is therefore a valid `SCIENTIFIC FAIL`. Keep
`pcg_pose_polish=no`. Do not weaken exact-pose, recovery, or FSC tolerances.
The isolated component is ready for scientific review: the remaining problem
is rotation stationarity/recovery against a reconstructed fixed map and noisy
FSC behavior, not the derivative, parallel loop, or execution harness.

**Status:** `ALL 16 ARMS COMPLETE; 64-THREAD EXECUTION PASS; SCIENTIFIC GATE
FAIL; FREQUENCY-LIMIT ARMS NON-DISTINCT`.

### Phase 11 — raw-data HPC beta-gal production A/B

**Disposition:** `DEFERRED — NOT CURRENT WORKFLOW EVIDENCE`. This historical
experiment ends with a direct `prg=reconstruct3D` call. Inspection of
`SIMPLE_data_testing/betagal.sh`, `apof.sh`, `motab.sh`, `proteasome.sh`, and
`trpm4.sh` found `abinitio3D` and `refine3D_auto` entry points but no direct
`reconstruct3D` entry point. The completed HPC checkpoint remains reusable if
the isolated component A/B is requested later, but this phase is removed from
the active acceptance path.

**Goal:** test whether an independently prepared experimental project shows the
same failure when the full SIMPLE workflow creates the starting checkpoint.
This test detects a bad frozen fixture or an incorrectly prepared starting
project, but it does not replace the truth-controlled Phase 10 diagnostic.

**Major file changes:**

- Required executable and validator preflight, fresh-run or resume selection,
  and guarded checkpoint handoffs:
  `.codex/hpc_betagal_prepare.sh:17`, `:19`, `:24`, `:51`, and `:66`.
- Existing-run resume mode selects a completed `*_refine3D*` project, records
  its checksum, and runs only the matched pose-polish validator:
  `.codex/hpc_betagal_prepare.sh:71` and `:78`.
- The conventional and hard-seeded `abinitio2D` comparison arms are preserved
  as comments because they do not test this 3-D end-polishing feature:
  `.codex/hpc_betagal_prepare.sh:125`, `:130`, and `:195`.
- Canonical beta-gal `abinitio2D -> model_cavgs_rejection ->
  abinitio3D_cavgs -> abinitio3D -> refine3D_auto` continuation from the
  untouched shared extraction project:
  `.codex/hpc_betagal_prepare.sh:205`.
- Matched shared/distributed `pcg_pose_polish=no|yes` validation from one final
  project, with logs and checksums kept under the experiment directory:
  `.codex/hpc_betagal_prepare.sh:249` and `:260`.

The existing raw-data preparation and 2-D hard-stream comparison remain
separate. The pose-validation source is the untouched shared extraction
project. Its scientific settings match `SIMPLE_data_testing/betagal.sh`:
`ncls=90`, `mskdiam=180`, `pgrp=d2`, and the standard five-stage 2-D-to-3-D
sequence without `inpl_cont` or `hard_stream`. Only `nthr` and `nparts` are
adjusted to fit the 48-CPU allocation. The final refinement project is frozen
once. The existing production validator then copies that checkpoint into
default, explicit-off, and enabled shared and distributed arms. Within each
matched pair, `pcg_pose_polish` is the only planned scientific difference.

Submit from the repository root on HPC:

```bash
sbatch .codex/hpc_betagal_prepare.sh
```

Do not repeat the expensive preparation when a completed project is available.
Downloaded experiment `Exp-20260817_112608` contains the completed checkpoint
`betagal-20260817_112608/10_refine3D_auto/betagal-20260817_112608.simple`.
Its pipeline log records `SIMPLE_REFINE3D NORMAL STOP`, cFAR `0.5626`, and a
4,633.6-second final refinement. The pose validator did not run because the old
checkpoint selector matched `_refine3D` but not `_refine3D_auto`. The selector
now accepts `*_refine3D*`. Resume only the production A/B with:

```bash
BETAGAL_REUSE_EXPERIMENT_ROOT=/home/hossainm7/Projects/Exp-20260817_112608 \
  sbatch .codex/hpc_betagal_prepare.sh
```

Resume mode does not rerun movie import, extraction, ab initio, or refinement.
It does not modify the frozen checkpoint. The production validator copies that
checkpoint into its own timestamped arms and compares default, explicit-off,
and enabled shared/distributed reconstruction paths. New validation output,
the checkpoint checksum, the resume log, and SLURM logs stay below the existing
`Exp-20260817_112608` directory.

All project outputs, the 3-D pipeline log, checkpoint checksum, validator log,
and nested `continuous_3D_pose_validation_*` directory remain under one
`/home/hossainm7/Projects/Exp-*` directory. Return that complete directory for
review. A normal SIMPLE exit is not sufficient. Review the nested
`analysis/summary.md`, pose errors, terminal accounting, FSC-area change, cFAR,
and shared/distributed agreement.

The apoferritin script is not changed in this phase. Its project uses
octahedral symmetry, while the current production analyzer supports only `c1`
and `d2`. Running it as `c1` or `d2` would make the pose-error comparison
scientifically invalid. Add octahedral symmetry support before using apoferritin
as a second dataset.

**Status:** `3-D PREPARATION COMPLETE; PRODUCTION A/B NOT RUN`. The downloaded
HPC evidence confirms the final refinement checkpoint. The previous validator
launch failure was a filename-pattern defect in the harness, not a scientific
failure. Resume support and the corrected selector are source-checked; the
resumed HPC execution and scientific review remain the user's validation gate.

## Decision hold for discussion with Hans

Phase 10C authorizes only the completed Phase 10D same-PCG-operator corrections.
Do not add the finite-box model $W$, change half-map ownership, alter LM
acceptance, or change the activation or workflow contract until Hans reviews
the truth matrix and the handoff. Keep `pcg_pose_polish=no` as the default. Do
not run Phase 11 as workflow evidence.

## Current next command

Phase 10F completed and failed the frozen scientific gate. Validate the
500-particle amortized batch with Oracle compilation, `case=pose_recovery`, the
mother suite, and a controlled production timing arm after synchronizing the
Oracle checkout. This is performance validation only and cannot reverse the
scientific result. Give the retained Phase 10F directory and
`continuous_3D_pose_end_polishing_handoff.md` to Hans and the refinement
researchers before changing the reference, rotation policy, frequency policy,
or acceptance rule. Frozen beta-gal and HPC A/B runs remain deferred.

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

The component implementation and production connection passed their original
focused, mother-suite, persistence, and shared/distributed gates, but the frozen
beta-gal A/B failed scientifically. Phases 10--10C separated the simulator
finite-box boundary, reconstructed-reference effects, and the KB-envelope
inverse crime. Phase 10D corrected the missing inverse-envelope application and
pose-only shell-range override. Oracle compilation, the envelope-free operator
gate, and the ten-case mother suite now pass. The clean repeated truth matrix
completed and failed exact-pose rotation stationarity, joint clean recovery,
and noisy FSC gates. It also found that the nominal frequency-limited arms did
not change the executed pose objective. The isolated investigation is now at
scientific handoff; the 500-particle batch is only a pending performance
validation. Frozen beta-gal and raw-data HPC direct-command A/B runs are no
longer workflow acceptance gates. Actual `abinitio3D` or `refine3D_auto`
integration, or a changed reference/frequency/acceptance policy, requires a new
or revised FINAL SPEC after Hans decides where and how often polishing belongs.
Stage 4 alternating reconstruction remains blocked.
