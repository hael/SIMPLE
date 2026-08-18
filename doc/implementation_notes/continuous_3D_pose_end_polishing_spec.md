# Continuous 3-D PCG post-reconstruction pose-polishing SPEC

**Status:** FINAL (operator-contract clarification, 2026-08-17)

**Frozen original proposal:** [continuous_3D_refinement_on_pcg_operator.md](continuous_3D_refinement_on_pcg_operator.md)
**PLAN:** [continuous_3D_pose_end_polishing_plan.md](continuous_3D_pose_end_polishing_plan.md)

## Decision summary

This FINAL SPEC is not awaiting interpretation or approval. It defines one
isolated, default-off `reconstruct3D` experiment. The implementation may be
committed while disabled because its activation, rollback, half-map ownership,
and operator contracts are explicit and its component tests pass.

The feature is **not scientifically accepted**. The completed truth matrix
failed exact-pose rotation stationarity, clean joint recovery, and noisy FSC
criteria. Therefore:

- keep `pcg_pose_polish=no` by default;
- do not recommend the option for production data;
- do not integrate it into `abinitio3D` or `refine3D_auto`;
- do not weaken the acceptance thresholds.

Hans and the refinement researchers are needed only if work continues under a
new or revised SPEC. That discussion must decide:

1. whether polishing is a terminal command, a repeated refinement stage, or
   should stop at the present research result;
2. what statistically valid fixed reference replaces or constrains the biased
   same-half reconstructed map;
3. whether reconstruction and polishing should share a new explicit frequency
   or rotation-regularization policy; and
4. what held-out or dataset-level signal may accept a particle update when the
   per-particle fixed-map objective decreases but truth pose or FSC worsens.

Everything below is the frozen technical contract for the implementation that
was tested. It is detail for verification, not an open design discussion.

## Request and outcome

The useful production unit is a joint local pose update, not a shift-only
update. Add one opt-in post-reconstruction pass to `reconstruct3D` that jointly
refines the three rotational and two image-shift coordinates against its fixed
PCG half-maps, then performs one final PCG reconstruction with the accepted
poses.

The public option is `pcg_pose_polish=yes`. Its default is `no`. It replaces the
uncommitted `pcg_shift_polish` prototype; there is no public shift-only mode.

## Activation contract

The polish is allowed only when all these conditions are true:

- `pcg_pose_polish=yes`;
- `rec_backend=pcg`;
- `combine_eo=no`;
- `prg=reconstruct3D`;
- the base reconstruction has produced independent even and odd PCG half-maps.

Invalid combinations must stop with a clear error. Internal workers and child
commands must not start a second polish. `refine3D` does not own a separate
pose-polishing implementation; a completed refinement project is polished by a
subsequent `reconstruct3D` invocation.

## Scientific contract

For particle $i$, first construct the fixed observed plane

$$
\bar y_i=Q_i(I_i;\eta_i),
$$

where $Q_i$ is the existing PCG particle-input path: noise normalization,
particle-edge taper, native FFT, packed-plane extraction, and whitening. The
statistics $\eta_i$ measured from $I_i$ are fixed data; they are not functions
of a trial pose or prediction. Minimize

$$
\Phi_i(R_i,t_i)=\frac{1}{2}\left\|
\bar y_i-\mathcal A_i(R_i,t_i)V_{h(i)}
\right\|_2^2 .
$$

Let $u_h$ denote the unconstrained solve coordinate and let $V_h=P u_h$ be the
already supported half-map written by PCG. Here $\mathcal A_i$ is the complete
executed observation operator, not a bare Fourier gather. Its required chain is

$$
\mathcal A_i(R_i,t_i)V_h=
T_i(t_i)\,G_P(R_i)\,E^{-1}V_h,
\qquad
T_i(t_i)=\frac{C_iS_i(t_i)}{\sqrt{\sigma_i^2}}.
$$

The factors have these fixed meanings:

- $P$ is the reconstruction support already represented in $V_h$. A soft
  support must not be applied a second time to the stored half-map;
- $E^{-1}$ is the PCG inverse Kaiser--Bessel gather envelope used before the
  forward gather when deapodization is enabled;
- $G_P(R_i)$ is the executed oversampled Fourier gather, with its actual fast
  interpolation value and derivative;
- $T_i$ uses the executed shift, CTF, `ctfflag`, whitening, sampling, precision,
  Fourier packing, and operation order from the PCG reconstruction operator;
- $Q_i$ is data-side preprocessing only. A mean, variance, edge profile, or
  amplitude statistic derived from the observed particle is measured once and
  held fixed while evaluating every trial objective.

This SPEC does **not** insert the simulator's padded inverse FFT, real-space
clip, or native FFT into $\mathcal A_i$. The current PCG reconstruction does not
apply that model-side operation. If evidence selects a finite-box model
$\mathcal B_i\mathcal A_i$, it must be introduced as a coordinated change to
the reconstruction forward operator, its adjoint and normal operator, and the
pose Jacobian. Adding $\mathcal B_i$ only to pose polishing is nonconforming.
A simulator-matched operator is not automatically the experimental-data
operator.

The pose Jacobian is the derivative of this complete chain. Since $P$,
$E^{-1}$, $C_i$, and $\sigma_i^2$ are independent of rotation, the rotation
columns differentiate $G_P(R_i)$ and then apply $T_i$; the shift columns
differentiate only $S_i(t_i)$. Finite differences of a gather that omits
$E^{-1}$, support, transfer, whitening, or the selected shell range do not
validate $\mathcal A_i$.

The local update has five coordinates,

$$
R_i\leftarrow R_i\exp([\omega_i]_\times),\qquad
t_i\leftarrow t_i+\delta t_i,qquad
q_i=(\omega_x,\omega_y,\omega_z,\delta t_x,\delta t_y).
$$

Requirements:

- use the executed fast Kaiser--Bessel value and derivative paths;
- verify that the fixed half-map is the supported PCG output, apply the PCG
  inverse envelope exactly once, and only then create the Fourier workspace;
- evaluate every predicted value and all five Jacobian columns with the exact
  PCG $T_iG_PE^{-1}$ chain from the already supported $V_h$; do not apply
  data-only normalization or taper to a trial prediction;
- use the PCG CTF flag, whitening, Fourier packing, phase, and precision
  conventions;
- use exactly the shell range used by the preceding PCG reconstruction. The
  reconstruction operator and pose workspace must obtain it from one shared
  setting; a pose-only shell restriction is not conforming;
- refine each particle only against its matching fixed half-map;
- keep state, half ownership, CTF metadata, observations, and search decisions
  unchanged;
- commit a pose only after a finite, fully recomputed objective improvement;
- restore the complete input pose for every rejected or unreliable result;
- scale radians and pixels explicitly, bound both step types, and report
  stencil switches and all terminal outcomes;
- treat weak or singular five-parameter systems as no-update results;
- keep the current local symmetry representative and evaluate pose error modulo
  the configured point group;
- anchor each local pass to its fixed half-map. Synthetic pose metrics must fix
  the global rigid-body gauge before comparison.

## Scope

- Generalize the validated shift component to a joint five-parameter Jacobian
  and per-particle Levenberg--Marquardt solve.
- Persist accepted rotations and shifts through the normal `ptcl3D` project
  path.
- Run one final PCG reconstruction after polishing.
- Keep shared and distributed scientific behavior equivalent.

## Out of scope

- Running a discrete or probabilistic orientation search.
- Alternating volume and pose updates; that remains Stage 4.
- CTF-parameter refinement or a new reconstruction backend.
- Platform-wide changes to the Kaiser--Bessel kernel.
- A claim that standalone PCG is superior to gridding.

## Acceptance criteria

1. Rotation and joint-pose derivatives agree with fixed-cell finite differences
   of the complete executed PCG operator $\mathcal A_i$, including support,
   inverse envelope, transfer, whitening, and shell restriction.
2. An independent forward generator passes a declared operator matrix that
   separates bare gather, inverse-envelope gather, finite-box gather, and their
   combination over realistic box/support margins. The matrix determines
   whether the independent generator is a valid test of the PCG model; it does
   not silently redefine that model. Any selected replacement production model
   and its adjoint must pass a numerical dot-product identity before entering
   PCG.
3. Known joint perturbations are recovered from multiple orientations; exact,
   weak, and singular cases remain unchanged.
4. Accepted objectives decrease monotonically, and all rotation and shift steps
   remain within their declared bounds.
5. Pose-error tests fix the global gauge and minimize over symmetry-equivalent
   rotations.
6. Even and odd particles use only their matching half-map.
7. Accepted rotations and shifts persist and are consumed by the final PCG
   reconstruction; disabled and invalid routes cannot change poses.
8. Shared and distributed reconstruction routes obey the same activation and
   ownership rules.
9. A frozen beta-gal A/B changes only `pcg_pose_polish`; angular changes remain
   local, terminal counts balance, and independent half-map FSC improves or is
   neutral under a predeclared tolerance.
10. Focused and mother tests pass on Oracle Linux.

## Constraints

- Preserve the original continuous-refinement implementation note unchanged.
- Preserve SIMPLE ownership boundaries and particle-cache policy.
- Oracle Linux is the compilation and runtime gate.
