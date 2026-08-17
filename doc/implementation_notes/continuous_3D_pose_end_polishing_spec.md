# Continuous 3-D PCG post-reconstruction pose-polishing SPEC

**Status:** FINAL

**PLAN:** [continuous_3D_pose_end_polishing_plan.md](continuous_3D_pose_end_polishing_plan.md)

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

For particle $i$, minimize the weighted fixed-half-map objective

$$
\Phi_i(R_i,t_i)=\frac{1}{2}\left\|
\frac{y_i}{\sqrt{\sigma_i^2}}-
\frac{C_i}{\sqrt{\sigma_i^2}}S(t_i)G(R_i)V_{h(i)}
\right\|_2^2 .
$$

The local update has five coordinates,

$$
R_i\leftarrow R_i\exp([\omega_i]_\times),\qquad
t_i\leftarrow t_i+\delta t_i,qquad
q_i=(\omega_x,\omega_y,\omega_z,\delta t_x,\delta t_y).
$$

Requirements:

- use the executed fast Kaiser--Bessel value and derivative paths;
- use the PCG CTF flag, whitening, Fourier packing, phase, and precision
  conventions;
- use the active reconstruction frequency limit, never a separate hand-set
  ladder;
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
   of the executed forward operator.
2. Known joint perturbations are recovered from multiple orientations; exact,
   weak, and singular cases remain unchanged.
3. Accepted objectives decrease monotonically, and all rotation and shift steps
   remain within their declared bounds.
4. Pose-error tests fix the global gauge and minimize over symmetry-equivalent
   rotations.
5. Even and odd particles use only their matching half-map.
6. Accepted rotations and shifts persist and are consumed by the final PCG
   reconstruction; disabled and invalid routes cannot change poses.
7. Shared and distributed reconstruction routes obey the same activation and
   ownership rules.
8. A frozen beta-gal A/B changes only `pcg_pose_polish`; angular changes remain
   local, terminal counts balance, and independent half-map FSC improves or is
   neutral under a predeclared tolerance.
9. Focused and mother tests pass on Oracle Linux.

## Constraints

- Preserve the original continuous-refinement implementation note unchanged.
- Preserve SIMPLE ownership boundaries and particle-cache policy.
- Oracle Linux is the compilation and runtime gate.
