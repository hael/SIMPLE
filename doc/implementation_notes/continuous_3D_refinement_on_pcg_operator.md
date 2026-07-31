# Continuous 3D refinement on the PCG operator: design note

## Status

**Planning only. No code is proposed for immediate implementation.**

This note expands the "Longer term" paragraph of section 10 of
`doc/policies/reconstruct3D_pcg_policy.md` from a statement of intent into an
implementable design. (It previously expanded section 12 of
`ctf_sigma_weighted_pcg_reconstruction.md`; that note was deleted and superseded
by the policy document in commit `5858d023f`, so its section numbering no longer
exists.) The policy note establishes the right high-level positions -- tangent-space
`SO(3)` increments rather than Euler angles, trust-region Gauss--Newton rather
than undamped Newton, continuous updates as a local refinement layer *after*
discrete search rather than a replacement for it, and half-set discipline.
Those positions are unchanged and are not re-argued here.

What this note adds is the part the policy note leaves open: the actual derivation of
the pose derivatives against *this* operator, the one new numerical primitive
they require, the per-particle solve structure, a cost model built on measured
numbers, a named hazard that the derivation exposes, and staged gates.

It also records an assessment of three attached sources (a technical brief on
SGD-based refinement, Zehni et al. 2019, and Brubaker et al. 2015), including
where this note deliberately departs from the brief's recommendation.

## 0. Position: what is already built

The validated `reconstructor_pcg` operator is the volume half of continuous
refinement, and it is the half that both reference papers implement least well.

| Component | Zehni et al. 2019 | Brubaker et al. 2015 | This codebase |
| --- | --- | --- | --- |
| Basis | Kaiser--Bessel window functions | sinc interpolation | Kaiser--Bessel (`kbinterpol`) |
| Volume update | ADMM on `½‖g − H(Θ)c‖² + λ‖∇c‖_{2,1}` | SAGD, line-searched Lipschitz constant | converged PCG on the normal equations |
| Conditioning | (handled inside ADMM) | first-order; ill-conditioning acknowledged | exact Hessian diagonal `rho`, shell-floored |
| Noise model | single global `σ²` | single global `σ` | per-particle, per-shell `σ²` |
| CTF | single PSF `h` | per-image `C_θ` | per-particle, `ctfflag`-aware |
| Support prior | TV | independent exponential, positivity | soft spherical support as `P H P` |

Zehni's density subproblem is, term for term, the system this operator solves.
Their `H(Θ)` is our `A`, and their ADMM inner solve is a linear system in
`H†H` -- which is `apply_normal` plus `solve`. TV-regularized reconstruction in
their sense is our existing solver plus a shrinkage step and a dual variable.

Section 4.2 of the attached brief argues at length that volume optimization is
ill-conditioned because low frequencies are oversampled and high frequencies
undersampled, and that one must precondition by the inverse diagonal of the
Hessian, estimated with Hutchinson's estimator. That argument is correct, and
it is an argument *for what is already built*: `build_precond` does not
estimate that diagonal stochastically, because `rho = Σ_i G_i† |T_i|²` **is**
that diagonal in closed form from the sampling geometry (exactly, up to the
banding the KB stencil introduces -- see section 5 of the policy note).

## 1. The one missing primitive

The forward model reads a Fourier plane sample at

```text
loc = padf * [h,k,0] . R
```

and interpolates the volume there with a 27-tap KB stencil (`gather_window`).
**The pose enters the model only through `loc`.** Therefore

\[
\frac{\partial}{\partial \text{pose}}\bigl(\text{plane value}\bigr)
= \nabla_{\text{loc}} \hat V \cdot \frac{\partial \text{loc}}{\partial \text{pose}} ,
\]

and `∇_loc V̂` comes from the *same* 27-tap stencil with the KB weights
replaced by their analytic derivatives.

That is the whole of it. `gather_window` is already the differentiable forward
model; it simply does not currently return the derivative. What is needed is
one new routine -- call it `dgather_window` -- returning the interpolated value
plus its three spatial derivatives.

No automatic differentiation framework is required, and none should be
introduced (section 2.3).

## 2. Assessment of the attached sources

### 2.1 Zehni et al. 2019 -- joint density and continuous angular refinement

Directly relevant, and the closest match to the architecture proposed here.
Their Algorithm 1 alternates an ADMM update of the density map with a gradient
step on continuous poses. Points that transfer:

- The alternating (block-coordinate) shape, which the policy note's section 11
  already specifies.
- The choice of KBWF for the basis, for the same stated reasons this codebase
  uses `kbinterpol`: compact support, isotropy, closed-form transforms.
- Empirical evidence that continuous refinement of poses initialized within
  roughly `σ_θ = 0.5 rad` recovers most particles.

Points that do **not** transfer, and matter:

- Their experiment uses `L = 500` synthetic projections of one molecule with a
  single global CTF and Gaussian noise. Per-particle CTF and per-particle,
  per-shell `σ²` -- the entire reason this operator exists -- are absent.
- **They report outliers, and attribute them correctly:** "the existence of
  outliers is explained by the non-convexity of the problem and getting stuck
  at local minima while optimizing for the 3D pose variables of some
  projection images." This is the central limitation of the whole approach and
  section 5.1 below is built around it.
- Their initial perturbation, `σ_θ = 0.5 rad ≈ 29°`, is a *coarse* regime. The
  basin of attraction narrows as resolution improves (section 5.1), so their
  demonstrated capture radius should not be read as applicable at high
  resolution.

### 2.2 Brubaker, Punjani and Fleet 2015 -- "Building proteins in a day"

Relevant as the origin of the SGD-for-cryo-EM claim, but note what the paper
actually does:

- It **marginalizes** poses over a quadrature grid on `SO(3)` and `R²`
  (their Eq. 4), with importance sampling to make the marginalization
  affordable. It does **not** perform continuous gradient-based pose
  refinement. Poses are latent variables integrated out, not parameters
  descended on.
- SAGD is applied to the **volume only**, with a line-searched Lipschitz
  constant and positivity enforced by truncation.
- Their prior is an independent exponential on voxel values.

So this paper is not evidence for direct gradient optimization of orientations.
It is evidence that (a) stochastic optimization over particles scales, and
(b) importance sampling makes marginalization affordable. Its genuinely
transferable idea is the importance-sampling trick, which is orthogonal to
continuous refinement and could equally accelerate the *existing* probabilistic
search machinery.

### 2.3 The SGD technical brief -- four disagreements

The brief recommends "Path A": adopt PyTorch/TensorFlow, represent orientations
as quaternions, build a differentiable forward model, and replace SIMPLE's
orientation search and reconstruction modules with a unified SGD loop. This
note recommends against all four elements, for the following reasons.

**(a) Do not introduce an autodiff framework.** The gradient required is a
three-term chain rule that is written out in full in section 3. Hand-writing it
is exact, cheaper than a taped autodiff graph, and keeps the build
Fortran-only. Introducing a Python autodiff runtime into a modern-Fortran HPC
application, to compute a derivative that fits on one page, is a large
architectural cost for negative benefit.

**(b) Do not replace the volume solve with preconditioned SGD.** The brief's
own section 4.2 establishes that first-order volume updates need Hessian-
diagonal preconditioning to converge at high resolution. This codebase already
has the exact diagonal and, more importantly, does not take gradient steps on
the volume at all -- it *solves* the normal equations, currently in 10-15 PCG
iterations at 1.53 s each. Replacing a converged preconditioned solve with
preconditioned first-order steps is a downgrade, not a modernization.

**(c) Do not replace the discrete orientation search.** Argued in section 5.1.
This agrees with the policy note's section 11, and disagrees with the
brief's section 7.1, which lists "the entire probabilistic orientation
assignment machinery ... must be replaced" as required work.

**(d) Prefer tangent-space increments to quaternions.** Quaternions use four
parameters for three degrees of freedom and require renormalization every step.
Parameterizing the *update* in the Lie algebra needs three unconstrained
parameters, has a trivial derivative at the origin, and -- decisively for this
codebase -- **requires no change to how orientations are stored**.
`prep_particles` already caches `rotmats(3,3,nptcls)`; only the update rule is
new. The policy note's section 11 already specifies this.

The brief is right that a resolution curriculum is essential -- it directly
mitigates the local-minimum failure Zehni et al. report -- but its prescription
("start at 60 A, then 40, then 20, coupled to a learning-rate schedule") is a
hand-set schedule, and this platform already has something better. Section 5.3
specifies mirroring `refine3D_auto` instead.

## 3. Forward model and derivatives

### 3.1 The model being differentiated

For particle `i` at plane index `(h,k)`, with `w_i = 1/sqrt(σ²_i(shell))`, the
whitened residual is

\[
r_i(h,k) = w_i\Bigl[\, C_i(h,k)\, S_i(h,k)\, \hat V\bigl(\mathrm{loc}_i(h,k)\bigr) - y_i(h,k) \Bigr],
\qquad
f = \tfrac12 \sum_i \sum_{h,k} \bigl| r_i(h,k) \bigr|^2 .
\]

`C_i` is the CTF (with `ctfflag` semantics as in `build_transfer`), `S_i` the
shift phase, and `V̂` the interpolated Fourier volume. This is the same
objective the volume solve minimizes; that is the point of doing it here rather
than in a separate engine.

### 3.2 Rotation derivatives

Parameterize the update in the tangent space at the current rotation:

\[
R \;\leftarrow\; R \exp\bigl([\omega]_\times\bigr), \qquad \omega \in \mathbb{R}^3 ,
\]

with `[ω]_× u = ω × u`. In this codebase's row-vector convention
`loc = padf · [h,k,0] · R`, so `loc(ω) = loc₀ · exp([ω]_×)` and, evaluating at
`ω = 0` with `G_a = [e_a]_×`,

\[
\frac{\partial\,\mathrm{loc}}{\partial \omega_a}
= \mathrm{loc}_0\, G_a
= -\,(e_a \times \mathrm{loc}_0)
= \mathrm{loc}_0 \times e_a ,
\qquad\text{so}\qquad
\frac{\partial\,\mathrm{loc}}{\partial \omega} = \bigl[\mathrm{loc}_0\bigr]_\times .
\]

**The Jacobian of the sample coordinate with respect to the rotation increment
is just the skew matrix of the current sample coordinate.** No matrix products,
no trigonometry, no chart bookkeeping -- `loc` is already computed in the
gather loop. Right-multiplication is the correct side: it applies the increment
in the volume frame, which is where the sampled coordinate lives.

Hence

\[
\frac{\partial r_i}{\partial \omega_a}
= w_i\, C_i S_i \,\bigl(\nabla \hat V \cdot (\mathrm{loc}_0 \times e_a)\bigr).
\]

Re-linearizing at the current `R` each outer iteration means gimbal lock and
chart singularities never arise.

### 3.3 In-plane rotation is not a separate degree of freedom

`omega` spans all of `so(3)`, so **the in-plane angle is already included and
requires no separate treatment**. Writing `delta_loc = loc x omega` and
splitting `omega` along the central-section normal `n = R(3,:)`:

- **`omega` parallel to `n`.** For `loc` lying in the section plane,
  `loc x omega` is in that same plane and perpendicular to `loc`. The section
  rotates within itself: this *is* the in-plane rotation.
- **`omega` perpendicular to `n`.** `loc x omega` acquires an out-of-plane
  component; the section tilts. These are the two projection-direction degrees
  of freedom.

The decomposition is useful for *interpreting* an update, for setting
per-component priors, and for diagnostics. It is **not** a computational split:
one `dgather_window` call and one 3x3 Jacobian `[loc]_x` produce all three
components together, at identical cost, with no special case for the in-plane
term.

This is the desired behaviour: continuous and uniform in every rotational
degree of freedom, on Cartesian central sections, with no polar resampling
anywhere in the refinement.

**Why the polar representation is not needed here.** SIMPLE's `pftc` machinery
exists because, in a *discrete correlation search*, an in-plane rotation is a
cheap index shift along the polar angle -- the right optimization when many
in-plane angles must be scored per particle. A gradient method never scores
many angles; it evaluates one derivative. The reason to separate in-plane from
projection direction therefore disappears entirely, and the Cartesian
central-section form the PCG operator already uses is the natural
representation. `pftc` is simply not involved in the polish stage, and no
on-the-fly polar extraction is required.

A consequence worth stating: the discrete search and the continuous polish will
run on **different Fourier representations** (polar for search, Cartesian for
polish). That is acceptable because the polish re-evaluates the objective in
its own representation and accepts or rejects on that basis, so no
cross-representation numerical agreement is assumed. It does mean the two
stages' objectives are not bitwise comparable, and stage-3 validation must
compare *before and after polish within the Cartesian model*, never
polish-objective against search-score.

### 3.4 Shift derivatives

The shift is the phase
`S_i(h,k) = exp(i·2π(h t_x + k t_y)/box)` (sign per `build_transfer`, which was
aligned to production's `gen_fplane4rec` convention). Therefore

\[
\frac{\partial S_i}{\partial t_x} = \frac{i\,2\pi h}{\text{box}}\, S_i ,
\qquad
\frac{\partial r_i}{\partial t_x} = w_i\, C_i \hat V(\mathrm{loc})\,\frac{i\,2\pi h}{\text{box}}\, S_i .
\]

Exact, analytic, and involving no interpolation at all. This is why shifts are
stage 1 (section 6): they exercise the entire outer-loop machinery on a
gradient that cannot be wrong.

### 3.5 The `dgather_window` primitive

`apod_mat_3d` and `apod_mat_3d_fast` build the 27-tap stencil as a **separable**
product of three 1-D weight vectors, each **normalized to unit sum**:

```text
kbw(i,j,k) = wx_hat(i) * wy_hat(j) * wz_hat(k),   wx_hat = wx / sum(wx)
```

Separability makes the derivative cheap: `∂kbw/∂loc_x = wx_hat'(i)·wy_hat(j)·wz_hat(k)`
and cyclically. Six 1-D vectors of length 3 are needed instead of three.

The 1-D normalization introduces a quotient-rule term that **must not be
dropped**. With `u_i = apod(base_x + i - 1)`, `base_x = win_lo_x − loc_x`,
`S = Σ_j u_j`:

\[
\frac{\partial \hat w_i}{\partial \mathrm{loc}_x}
= -\frac{1}{S}\left( a'_i - \hat w_i \sum_j a'_j \right),
\qquad a'_i \equiv \mathrm{apod}'(\text{base}_x + i - 1).
\]

`apod(x) = (1/W)·I₀(β√(1−(2x/W)²))` for `|x| ≤ Whalf`, so

\[
\mathrm{apod}'(x) = -\frac{1}{W}\, I_1\!\bigl(\beta\sqrt{u}\bigr)\, \frac{\beta}{\sqrt{u}} \cdot \frac{4x}{W^2},
\qquad u = 1 - (2x/W)^2 .
\]

`bessi1` already exists in `simple_kbinterpol` (module-private), which is a
strong argument for the derivative routine living **in that module** as a
sibling of `apod_mat_3d_fast`, rather than in the reconstruction module.

Cost: the weight computation roughly doubles; the 27-tap accumulation goes from
one sum to four. Call it ~3x a plain gather.

### 3.6 Per-particle Gauss--Newton

With 3 rotation parameters (or 5 including shifts) per particle, the problem is
**separable per particle and tiny**. Assemble, per particle, over its plane:

\[
(J^H J)_{ab} = \Re \sum_{h,k} \overline{\partial_a r}\; \partial_b r ,
\qquad
g_a = \Re \sum_{h,k} \overline{r}\; \partial_a r .
\]

Solve the damped 5x5 system `(J^H J + μ diag) δ = −g` (Levenberg--Marquardt),
apply `R ← R exp([δ_ω]_×)` and `t ← t + δ_t`, and accept or reject on a
recomputed objective. Two or three such steps converge.

This is a deliberate departure from the brief's SGD-with-schedules
recommendation. Learning rates, momentum and annealing are the right tools when
parameters are numerous, coupled, and the Hessian cannot be formed. Here there
are five parameters per particle, they are uncoupled across particles, and the
Gauss--Newton Hessian is a 5x5 matrix obtained for free from derivatives
already computed. Importing schedule-tuning machinery would be solving a
problem this formulation does not have.

This satisfies the standing requirement that the implementation expose the
residual, Jacobian-vector product and adjoint-Jacobian-vector product as tested
operator methods, and never form a dense particle-parameter Hessian. A
per-particle 5x5 is not one: the full matrix remains block-diagonal and is
never assembled.

### 3.7 Consistency requirements

The pose refinement must use **the same forward model as the reconstruction**,
or it will optimize alignment under one convention and reconstruct under
another -- the exact failure the policy note's section 11 warns about.
Concretely:

- **Deapodization.** The reconstruction solves with `A = Ã·E⁻¹`. The predicted
  data for particle `i` is therefore `Ã_i(E⁻¹x)`. Since `E⁻¹` does not depend
  on pose, this reduces to deapodizing once and gathering from the result --
  which is what `apply_normal_matrixfree` already does via `set_volume`.
- **`ctfflag`.** `CTFFLAG_NO` and `CTFFLAG_FLIP` must be honoured exactly as in
  `build_transfer` and `absT2_plane`. Getting this wrong is invisible on
  synthetic fixtures (which default to `CTFFLAG_YES`) and destroys real-data
  results; this exact bug occurred once in this codebase and was fixed during
  the milestone-3 pass.
- **Support mask.** If the volume was solved under `P H P`, the predicted data
  must be generated from the masked volume.
- **Oversampling.** All coordinates are on the `padf`-times lattice.
- **Resolution limit.** The pose polish must mask on the same `sqlp` as the
  volume solve in the same outer iteration, driven by `params%kfromto(2)`
  (section 5.3). Refining poses against shells the volume solve excluded
  fits noise directly.

## 4. Named hazard: the interpolant is not C¹

Deriving the gradient exposed a property of the existing interpolation that is
harmless for one-shot reconstruction and is **not** harmless for gradient-based
pose refinement.

`apod(x)` is `(1/W)·I₀(β√(1−(2x/W)²))` for `|x| ≤ Whalf` and identically zero
beyond. At `x = Whalf` the argument of `I₀` is zero, so `apod(Whalf) = 1/W`
-- **not** zero. For the configured `KBWINSZ = 1.5, KBALPHA = 2`, `apod(0) ≈ 35.3`
and `apod(1.5) ≈ 0.333`, so the kernel steps to zero from about **1% of its
peak value**.

The stencil corner is `win_lo = nint(loc) − iwinsz`, which is piecewise
constant in `loc`. As `loc` crosses a half-integer the window shifts by one
tap: the tap sitting at `−Whalf` is dropped and a new one appears at `+Whalf`.
The interpolated value therefore has a small jump discontinuity at every
stencil switch, of relative order 1%, and the derivative has a corresponding
jump.

Consequences, and what to do:

- The objective is **piecewise smooth**, not smooth. Naive line searches and
  strict Wolfe conditions will misbehave.
- This is a positive argument for the trust-region accept/reject step that
  the policy note's section 11 already mandates: acceptance is decided on a
  recomputed
  objective, which is robust to a non-differentiable point being crossed.
- **Oversampling is a real mitigant.** On the `padf`-times lattice `V̂` is
  band-limited to half the padded Nyquist and hence much smoother between
  samples, so the value mismatch across a switch is far smaller than the raw
  1% figure suggests. Quantifying this is a stage-2 measurement, not a guess.
- A cleaner long-term fix is a KB kernel that tapers to zero at `|x| = Whalf`,
  but that changes the interpolation used by the whole platform and is
  explicitly **out of scope** for this work.

Stage 2's gate must therefore include measuring the actual size of the jump
across a stencil switch on the padded lattice, not only a finite-difference
check away from one.

## 5. Optimization shape and cost

### 5.1 Local refinement, not global search

The distinction that matters is not gradient versus expectation-maximization;
it is **local versus global**. Neither attached paper replaces global search
with gradients: Zehni et al. refine locally from a perturbed initialization and
report local-minimum outliers; Brubaker et al. marginalize over a quadrature
grid.

The basin of attraction narrows as resolution improves. A rotation error `δ`
displaces a feature at radius `R` pixels by `R·δ` pixels. Requiring that
displacement to stay below one resolution element `d_pix = d_Å / smpd` gives

\[
\delta \lesssim \frac{d_{\text{pix}}}{R} .
\]

For the bgal test case (`box 256`, `mskdiam 180`, so `R ≈ 90 px`) at 3 Å with
`smpd ≈ 1 Å/px`, this is `δ ≲ 0.033 rad ≈ 1.9°`, and one would want a comfortable
margin below that. Covering the full sphere at spacing `δ` needs roughly
`4π/δ² ≈ 1.1 × 10⁴` directions at 1.9°, and about `1.6 × 10⁵` at 0.5°.
(Order-of-magnitude only, and it concerns the cost of the *discrete* search:
SIMPLE searches the asymmetric unit and scores in-plane rotations separately in
the polar representation, so the operative `nspace` is smaller. No such
separation exists in the continuous formulation -- see section 3.3.)

**This is the actual argument for continuous refinement, and the brief does not
make it crisply:** the cost of a grid search grows as `nspace × nptcls` with
`nspace ∝ 1/δ²`, whereas the cost of continuous refinement is *independent of
the target angular precision*. Continuous refinement is not a better search --
it is the removal of discretization error that search cannot remove without
quadratic cost.

The design consequence is that this should be **a polishing stage appended to
the existing pipeline**, not a replacement for it. Global search puts each
particle in the right basin; gradients remove the residual grid error. This is
also what production packages do (RELION local searches, cryoSPARC local and
non-uniform refinement).

### 5.2 Cost model

Measured on 5000 particles, `box 256`, 10 threads:

| Quantity | Measured |
| --- | --- |
| fused gather+weight+scatter particle pass | 9.93 s |
| kernel operator, per PCG iteration | 1.53 s |
| matrix-free operator, per PCG iteration | 11.5 s |
| setup (I/O, RHS, one particle pass, kernel build) | 26.9 s |
| full volume solve, 10 iterations | ~53 s total |

Estimates that follow, flagged as estimates:

- A gather-only pass is roughly half the fused pass, so ~5 s; a derivative
  gather at ~3x is **~15 s per outer iteration** for pose refinement over all
  particles. The per-particle 5x5 solves are negligible.
- **The kernel must be rebuilt whenever poses change.** `Khat` and `rho` both
  depend on the rotations, so an alternating scheme pays the full setup cost
  every outer iteration. Budget roughly `27 s (rebuild) + 25 s (solve) + 15 s
  (pose polish) ≈ 65 s` per outer iteration, not the 25 s that the kernel's
  cheap iterations might suggest in isolation.

This weakens, but does not remove, the kernel's advantage in an alternating
scheme, and it is the strongest argument for eventually restricting pose
refinement to minibatches (section 5.4) or for accepting fewer, better outer
iterations.

### 5.3 Resolution limit: mirror `refine3D_auto`, do not invent a schedule

The resolution curriculum is **not** a hand-set schedule and must not become
one. `refine3D_auto` already drives it as a closed loop through the nonuniform
filter, and this design mirrors that working approach rather than replacing it.

The mechanism, as implemented:

- `refine3D_auto` sets `filt_mode = 'nonuniform'`, `nu_refine = 'yes'`,
  `lplim_crit = 0.143` and `incrreslim = 'no'`
  (`simple_commanders_refine3D.f90`).
- The starting limit is seeded from the *starting volume's own FSC* at the
  0.143 criterion, in `seed_refine3D_auto_nonuniform_lpset`, which sets
  `params%lp` and `params%kfromto(2)`.
- Every subsequent iteration re-derives the limit in `set_bp_range3D`
  (`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`):
  `try_nu_refine_project_alignment_lp` returns the NU filter's own matching
  low-pass when available, otherwise `get_find_at_crit(fsc, lplim_crit)` takes
  it from the current FSC; `lpstop` clamps it; it is then clamped to the
  cropped Nyquist and written back to the project with
  `spproj_field%set_all2single('lp', params%lp)` for the next iteration to read.

So the search limit tracks the resolution the reconstruction has *actually
achieved*, and frequency marching is emergent rather than prescribed. That is
strictly better than any fixed ladder, and it removes what would otherwise be
this design's most arbitrary tuning surface. It also means the NU filter's
local-resolution estimate -- not a global scalar -- is what ultimately governs
the limit, which matters for anisotropically resolved maps.

Implications for this work:

- The operator's resolution limit is currently fixed at native Nyquist
  (`self%sqlp = R*R` in `new`). It needs a setter so the outer loop can drive it
  from `params%kfromto(2)`, exactly as the matcher does.
- `absT2_plane` masks on `sqlp`, so `rho` and `Khat` both depend on the limit.
  **Changing the limit forces a preconditioner and kernel rebuild** -- the same
  rebuild already forced by moved poses (section 5.2), so the two costs
  coincide rather than compound.
- The pose polish must use the *same* limit as the volume solve in the same
  outer iteration. Refining poses against shells the volume solve excluded
  would fit noise directly.
- Nothing here should be reimplemented. The outer loop should consume
  `params%kfromto(2)` and the project's `lp` field, and let the existing NU
  machinery own the policy.

### 5.4 Why not SGD, yet

The full volume solve is ~53 s for 5000 particles. For typical datasets a batch
alternating scheme may simply be fast enough, and the policy note's section 11
already sets out that minibatch stochastic optimisation is a later option with
PCG remaining the reference volume solver. Minibatching adds step-size schedules, variance
control, and stratification requirements (half sets, orientation coverage,
optics groups) -- real complexity that should be paid for only once a measured
need exists. Defer to stage 5 and revisit with numbers.

## 6. Staged implementation with gates

Each stage must pass its gate before the next begins. This mirrors the staging
that made the reconstruction operator trustworthy.

**Stage 0 -- even/odd split and FSC in `reconstruct3D_pcg`.**
This comes first and is not optional. At present there is *no quantitative way
to tell whether any change helps* -- every judgement so far has been visual
inspection in Chimera. Before adding thousands of free parameters to the fit, a
gold-standard resolution curve must exist. This is also a prerequisite for the
real-space priors sketched in the policy note's section 11 and for the half-set
discipline in section 7.
*Gate:* two independent half-set volumes and an FSC curve, reproducing a
sensible resolution on a known dataset.

**Stage 1 -- shifts only.**
Analytic gradient (section 3.4), no new interpolation primitive, 2 parameters
per particle. Validates the entire outer-loop scaffold -- residual assembly,
Gauss--Newton solve, trust-region accept/reject, convergence monitoring -- on a
gradient that cannot be wrong.
*Gate:* synthetic recovery of known injected shifts; FSC improves or is neutral
on real data.

**Stage 2 -- `dgather_window`.**
The one genuinely new piece of numerics (section 3.5), placed in
`simple_kbinterpol` beside `apod_mat_3d_fast`.
*Gate:* finite-difference agreement against `forward_plane`, in the same spirit
as the adjoint dot-product test that gated the operator -- that test is what
made the operator trustworthy and this is its counterpart. **Plus** measurement
of the stencil-switch discontinuity from section 4.

**Stage 3 -- rotation refinement.**
Tangent-space increments covering all three rotational degrees of freedom
including in-plane (section 3.3), per-particle Levenberg--Marquardt, local
only, initialized from the existing discrete assignment. Resolution limit taken
from `params%kfromto(2)` as driven by the NU filter (section 5.3) from the
outset -- never a hand-set ladder.
*Gate:* synthetic recovery of known injected rotations; on real data, angular
changes are small (consistent with grid spacing, not wholesale reassignment)
and half-set FSC improves.

**Stage 4 -- alternate.**
Volume solve (existing) alternating with pose polish. This is Zehni's Algorithm
1 with a stronger solver on both halves.
*Gate:* outer-loop convergence is monotone in the recomputed objective; FSC
improves over the fixed-pose baseline; no drift in particle rejection.

**Stage 5 -- minibatching, only if a measured need exists.**
See section 5.3 above and the policy note's section 11.

## 7. Validation and identifiability

Section 12.3 lists the required safeguards and they stand. Two points deserve
sharpening in light of the design above.

**Overfitting arithmetic.** Stage 3 introduces 3 parameters per particle, stage
4 with shifts introduces 5. At 5000 particles that is 25,000 free parameters
fitted against the data, alongside a 16.7 M-voxel volume. This is precisely the
regime where a resolution estimate can measure one's own overfitting.

**Half-set discipline is a design constraint on stage 0, not an afterthought.**
Poses refined against a half-map must be applied only to that half. If a
particle's orientation is refined against a volume built partly from that same
particle, the gold-standard FSC is no longer gold-standard, and the refinement
will appear to improve resolution precisely because it is fitting noise. The
even/odd split must therefore be in place *before* stage 1, and the outer loop
must carry two independent (volume, pose-set) pairs throughout.

**Numerical precision.** Pose gradients near convergence involve differences of
nearly equal quantities. `solve` already accumulates its dot products in double
precision; the derivative gather and the per-particle normal equations will
likely need the same treatment. Worth measuring at stage 2 rather than assuming
either way.

## 8. Deliberate non-goals

- No autodiff framework, no Python/C++ runtime dependency.
- No replacement of the discrete/probabilistic orientation search.
- No first-order (SGD/SAGD/Adam) update of the volume; PCG remains the volume
  solver and the reference.
- No change to orientation storage; `rotmats` stays as it is.
- No hand-set resolution schedule, and no reimplementation of low-pass policy.
  The limit comes from `params%kfromto(2)` as the NU filter drives it
  (section 5.3).
- No polar (`pftc`) representation in the polish stage, and no separate
  treatment of the in-plane angle (section 3.3).
- No change to the platform-wide KB kernel (section 4).
- No per-particle CTF refinement at this stage; the policy note stages that
  separately and it should stay staged separately.

## 9. Open questions

1. **Does the kernelized operator remain valid as poses move?** It was
   validated against the matrix-free reference at fixed poses. Its
   shift-invariance approximation may interact with pose updates in ways not
   yet probed. Stage 4 should compare an outer loop run with `pcgop=kernel`
   against one with `pcgop=matrixfree`.
2. **How large is the stencil-switch discontinuity on the padded lattice, in
   practice?** Section 4 argues oversampling suppresses it; that is an
   argument, not a measurement.
3. **Is per-particle Levenberg--Marquardt or a shared trust region better?**
   Per-particle damping is simpler; a shared radius couples the particles and
   may be more stable early in frequency marching.
4. **Does the NU-driven limit remain well behaved once poses move?**
   Section 5.3 settles the *policy* question -- the limit is taken from the NU
   filter exactly as `refine3D_auto` does it, not invented here. What is not
   settled is the feedback: in `refine3D_auto` the limit responds to a
   reconstruction whose poses came from a discrete search, whereas here poses
   and limit would co-adapt. A limit that advances because poses overfit, which
   then licenses further overfitting, is a plausible failure mode. The half-set
   discipline in section 7 is the primary guard; whether an additional cap on
   per-iteration limit advance is needed is an open empirical question for
   stage 4.
5. **Should importance sampling (Brubaker et al. section 3.3) be applied to the
   existing probabilistic search** instead of, or before, any of this? It is
   orthogonal to continuous refinement and might deliver a larger speedup for
   less architectural risk.

## References

- Zehni, Donati, Soubies, Zhao, Do, Unser. *Joint density map and continuous
  angular refinement in cryo-electron microscopy.* IS&T Electronic Imaging,
  Computational Imaging XVII, 133-1 (2019).
- Brubaker, Punjani, Fleet. *Building proteins in a day: efficient 3D molecular
  structure estimation with electron cryomicroscopy.* arXiv:1504.03573 (2015);
  IEEE TPAMI 39(4) (2017).
- Punjani, Rubinstein, Fleet, Brubaker. *cryoSPARC.* Nature Methods 14(3),
  290-296 (2017).
- Scheres. *RELION.* J. Struct. Biol. 180(3), 519-530 (2012).
- `doc/policies/reconstruct3D_pcg_policy.md` -- contract of the code that runs
  today; especially section 5 (preconditioner and kernelized operator) and
  section 10 (deferred and future work). Supersedes the deleted
  `ctf_sigma_weighted_pcg_reconstruction.md`.
- `doc/implementation_notes/pcg_reconstruction_as_gridding_replacement.md` --
  workflow integration; its section 5.3 establishes that the PCG *solve* is
  data-free, which does NOT hold for the pose objective specified here.
