# Proposal: complementary molecular and solvent priors for PCG 3D refinement

## Status

**Research proposal only. No implementation is included in this note.**

The proposed features are opt-in and disabled by default. The first adds a
convex real-space solvent-flatness penalty to the existing PCG reconstruction
objective. The second adds a molecular spectral prior based on Wilson
statistics. Both use the state-specific nonuniform (NU) molecular envelope as
a graded partition. Neither hard-zeros the solvent region.

The main recommendation is:

> Keep the existing broad spherical support constraint, and add a weighted
> solvent-variance precision term to the PCG normal operator used for the one
> even reconstruction and the one odd reconstruction. A Wilson molecular prior
> can be added in the same solves if it is represented by a fixed quadratic
> precision operator. Do not insert literal LocScale amplitude rescaling into
> conjugate gradients.

This preserves the current one-solve-per-half workflow. The prior changes the
half maps, and their FSC can therefore change when the feature is enabled. It
does not change particle accumulation or the data likelihood.

## 1. Decision requested

Approve a staged experiment with the following fixed design:

1. Add a dimensionless control named `pcg_solvent_lambda_rel` with default
   value `0.0`.
2. Interpret a positive value as the strength of a graded solvent-flatness
   prior relative to the PCG data scale.
3. Use the previous available `nu_envmask3D_stateNN.mrc` as a molecular
   envelope. Use its complement as solvent confidence.
4. Hold the envelope fixed during each PCG solve and use the same envelope for
   both half sets.
5. Add the solvent precision before the normal even and odd solves. Do not add
   another solve, warm start, raw-accumulator reload, or `_unfil` product for
   this feature.
6. Calculate the ordinary FSC from the resulting even and odd maps, then run
   the existing NU filtering and envelope-generation steps.
7. Treat the envelope generated in the current iteration as input to a later
   iteration, never as mutable state inside the current PCG solve.
8. Skip the envelope-dependent priors for that iteration if the NU envelope is
   missing or invalid. Do not silently substitute a density mask or a
   spherical mask.
9. Design the PCG prior interface so a fixed Wilson molecular precision can be
   added independently with `pcg_wilson_lambda_rel`, also defaulting to `0.0`.
10. Test four separate cases: neither prior, solvent only, Wilson only, and
    both. Do not infer the value of the combination from either prior alone.

## 2. Motivation and literature position

LocScale and solvent flattening solve related but different problems.
LocScale-1.0 performs local radial Fourier-amplitude scaling against a model
map [2]. The physics-informed route in LocScale-2.0 constructs a pseudoatomic
or hybrid reference and scales local Fourier amplitudes while preserving the
observed local phases [1]. Its scattering expectations are connected to
Wilson statistics [3-6]. LocScale-2.0 lists solvent flattening as a possible
*future* real-space prior; it does not include solvent flatness in its current
physics-informed reconstruction objective [1].

Wilson's 1942 paper explains why the averaged high-frequency intensity of a
random collection of atoms approaches the sum of atomic intensities [3]. His
1949 paper derives probability distributions for crystallographic
intensities [4]. Singer gives a modern derivation and states the regime in
which Wilson statistics is valid [5]. Gilles and Singer use the resulting
molecular mean and covariance to define a Bayesian molecular prior [6]. These
papers are important for a future molecular-spectrum or covariance prior, but
they are not the original source of solvent flattening.

Solvent flattening starts with Wang's iterative use of a molecular boundary
and the expectation of nearly constant density in the solvent [7]. Cowtan and
Main combine real-space and reciprocal-space constraints [8]. Terwilliger
expresses solvent flattening as an explicit likelihood problem [9]. The later
cryo-EM density-modification work uses expectations for both solvent and
macromolecular regions [10]. This is the direct conceptual line for the
proposed PCG penalty.

The proposal therefore keeps three ideas separate:

- **Wilson prior:** a molecular mean/covariance model.
- **LocScale operation:** local, phase-preserving amplitude rescaling against
  an expected or reference spectrum.
- **Solvent prior:** expected low real-space variance in the solvent.

Wilson and solvent precisions can be combined in one linear PCG system.
LocScale can contribute a fixed local covariance profile, but literal amplitude
matching needs a nonlinear outer iteration. All components must first be tested
individually before the joint result can be interpreted.

## 3. Goals and non-goals

### 3.1 Goals

- Replace hard solvent truncation with a graded penalty.
- Preserve the PCG operator's self-adjoint, positive-semidefinite form.
- Reuse the existing lagged NU envelope in the primary half-map solves.
- Preserve one reconstruction solve for each half.
- Preserve kernel and matrix-free operator parity.
- Add only linear-time real-space work to each PCG matrix-vector product.
- Make missing-mask and weak-mask behavior explicit and diagnosable.

### 3.2 Non-goals for the first experiment

- Do not remove the broad spherical support constraint.
- Do not implement full crystallographic density modification.
- Do not add histogram matching, positivity, solvent flipping, or a learned
  denoiser.
- Do not reproduce LocScale's windowed amplitude-scaling pipeline inside the
  reconstruction operator without a fixed quadratic surrogate.
- Do not conflate a Wilson covariance prior with LocScale postprocessing.
- Do not update the NU envelope during a PCG solve.
- Do not add a replay solve or a solvent-specific `_unfil` map pair.
- Do not change particle weights, raw `(B,D)` statistics, or trailing-chain
  semantics.

## 4. Proposed mathematical model

### 4.1 Existing PCG estimator

Let `H_data` be the data-only normal operator, `b` the adjoint-data right-hand
side, and `lambda_0` the small PCG base Tikhonov term. Let `P` be the existing
diagonal soft spherical support operator. The solver represents the map as

\[
x = P u.
\]

Before adding the solvent prior, the per-half solve considered here has the
form

\[
P\left(H_{\mathrm{data}} + \lambda_0 I\right)P u = P b.
\]

The current-iteration FSC does not exist until both half solves finish. It is
an output of this reconstruction sequence, not an input to the solvent prior.

The concrete matrix-free and kernel operators remain data/prior operators.
The shared `apply_normal` wrapper applies the outer `P H P` contract.

### 4.2 Graded solvent confidence

Let `m_v` be the NU molecular envelope at voxel `v`, with

\[
0 \le m_v \le 1,
\]

where `m_v = 1` means confident molecular region and `m_v = 0` means confident
solvent. Let `p_v` be the existing broad support value. Define

\[
w_v = p_v(1-m_v), \qquad s_v = w_v^2, \qquad
S = \sum_v s_v.
\]

Including `p_v` in the weight prevents the many fixed zero voxels outside the
broad reconstruction support from dominating the solvent mean. It also grades
the prior smoothly at both the NU-envelope boundary and the spherical-support
boundary.

If the NU machinery supplies solvent confidence directly, the adapter can use
that confidence in place of `1-m_v`. The stored artifact must nevertheless
have one documented convention.

### 4.3 Weighted solvent mean and penalty

For a trial map `x`, define the weighted solvent mean

\[
\mu_s(x) = \frac{\sum_v s_v x_v}{S}.
\]

The proposed penalty is

\[
R_s(x) = \frac{\lambda_s}{2}
\sum_v s_v\left[x_v-\mu_s(x)\right]^2.
\]

This is the graded form of solvent flattening. Voxels with high solvent
confidence receive the full variance penalty. Voxels in the soft transition
receive less. Confident molecular voxels receive none.

The user's suggested operator can be written precisely as

\[
L = W C_s, \qquad W=\operatorname{diag}(w), \qquad
C_s = I - \mathbf{1}\frac{s^T}{S}.
\]

Here `C_s x = x - mu_s(x) 1`. The normal operator must add
`L^T L`, not `L` itself. `L` is generally not symmetric because weighted
centering is not an orthogonal projection in the ordinary Euclidean inner
product. Adding `L` directly would invalidate the symmetry assumption used by
conjugate gradients.

The required precision operator is

\[
Q_s = L^T L
    = D_s - \frac{s s^T}{S},
\qquad D_s=\operatorname{diag}(s),
\]

and its matrix-free application is

\[
Q_s x = s \odot \left[x-\mu_s(x)\right].
\]

No matrix is stored. One weighted sum and two volume passes are sufficient.

### 4.4 Proposed per-half normal equation

The single even or odd solve becomes

\[
P\left(H_{\mathrm{data}} + \lambda_0 I +
\lambda_s Q_s\right)P u = P b.
\]

The effective strength is tied to the existing deterministic data scale:

\[
\lambda_s =
\texttt{pcg_solvent_lambda_rel}\;s_{\mathrm{data}}(D).
\]

This makes the new input dimensionless and preserves its meaning when data are
partitioned, reduced, fractionally updated, or uniformly reweighted.

### 4.5 Required algebraic properties

For every real map `x`,

\[
x^T Q_s x = \sum_v s_v[x_v-\mu_s(x)]^2 \ge 0.
\]

Therefore `Q_s` is symmetric positive semidefinite. It has the intended null
space:

\[
Q_s(c\mathbf{1})=0.
\]

It penalizes solvent variation, not the solvent offset. This is important
because the absolute map offset is normally not determined by cryo-EM data.

The broad support `P` should remain in the first implementation. The solvent
operator has an unpenalized constant mode, and the PCG interpolation/deapodize
path is least reliable near the padded-box corners. Removing `P` would add many
weakly constrained unknowns and change conditioning at the same time as the
new prior is evaluated.

## 5. Can the molecular and solvent priors be applied together?

### 5.1 Short answer

Yes. The combination is scientifically meaningful because the priors encode
different information in complementary regions:

- the solvent prior says that density outside the molecular envelope should
  have low variance around an unknown constant;
- a Wilson prior says how molecular scattering potential is distributed over
  spatial frequencies and, in its fuller form, how nearby Fourier
  coefficients are correlated [5,6].

The combination is mathematically compatible with PCG if both priors are held
fixed during a half solve and each is expressed as a quadratic positive-
semidefinite precision term. The sum of such terms remains symmetric positive
semidefinite.

This is best viewed as a product-of-experts model. The likelihood, molecular
prior, and solvent prior contribute additive negative log probabilities. It
does not require the two priors to be literally independent. Their separate
strengths state how much confidence is assigned to each expert.

### 5.2 A quadratic Wilson molecular prior

Let

\[
M=\operatorname{diag}(m)
\]

be the soft molecular selector from the NU envelope. Let `mu_W` and `C_W` be a
fixed Wilson-derived mean and covariance for the molecular scattering
potential, and let

\[
Q_W=C_W^{\dagger}
\]

be its positive-semidefinite precision on the supported subspace. A Gaussian
Wilson prior can be written as

\[
R_W(x)=\frac{\lambda_W}{2}
\left(Mx-\mu_W\right)^T Q_W\left(Mx-\mu_W\right).
\]

Singer derives the non-diagonal covariance of Fourier coefficients from the
random bag-of-atoms model [5]. Gilles and Singer use the corresponding mean and
covariance to construct a molecular prior for linear inverse problems [6].
This is the form that transfers naturally to PCG. It is a covariance model,
not a sharpening operation.

For a first implementation, a shell-diagonal approximation is the lowest-risk
choice:

\[
Q_W = F^*\operatorname{diag}(q_W(k))F,
\]

where `q_W(k)` is the inverse expected molecular variance, regularized away
from zero and normalized over the active reconstruction band. This corresponds
to the tractable radial covariance approximation used in conventional Bayesian
refinement. A later implementation can test the short-range Fourier
correlations described by Singer and Gilles [5,6].

A zero-mean choice, `mu_W=0`, adds only precision and suppresses coefficients
according to their expected variance. It does **not** force the reconstruction
to have a particular nonzero amplitude spectrum. A nonzero Wilson mean adds a
fixed prior contribution to the right-hand side and therefore needs stronger
provenance and bias controls.

LocScale's expected local power profiles can inform the covariance without
performing amplitude replacement. For fixed window extractors `R_j` and fixed
inverse-variance spectra `q_{W,j}`, a local molecular precision is

\[
Q_{W,\mathrm{local}}=
\sum_j M R_j^T F^*\operatorname{diag}(q_{W,j})F R_j M.
\]

Every summand is positive semidefinite, so the sum is PCG-compatible while the
profiles remain fixed. This is a LocScale-inspired **local covariance prior**,
not the LocScale scaling algorithm. It is more expensive than the global
shell-diagonal approximation and belongs after the global form is validated.

### 5.3 Combined PCG objective

With `x=Pu`, the combined objective is

\[
\begin{aligned}
J(u) ={}& J_{\mathrm{data}}(Pu)
 + \frac{\lambda_0}{2}\|Pu\|^2 \\
&+ \frac{\lambda_s}{2}(Pu)^TQ_s(Pu) \\
&+ \frac{\lambda_W}{2}
  \left(MPu-\mu_W\right)^TQ_W
  \left(MPu-\mu_W\right).
\end{aligned}
\]

The normal equation is

\[
P\left[H_{\mathrm{data}}+\lambda_0I+\lambda_sQ_s+
\lambda_W M Q_W M\right]Pu
=P\left[b+\lambda_W M Q_W\mu_W\right].
\]

For `mu_W=0`, the right-hand side remains the data-only `b`. Both new terms are
applied in the same even or odd solve. No replay is introduced.

The envelope provides a soft division of responsibility. `M` localizes the
molecular prior, while `w=p(1-m)` localizes solvent flatness. Both terms act in
the transition zone, but with reduced weight. This overlap is intentional: a
hard partition would introduce boundary ringing and make the result sensitive
to one mask contour.

The effective strengths should be independent:

\[
\lambda_s=\alpha_s s_{\mathrm{data}}(D),\qquad
\lambda_W=\alpha_W s_{\mathrm{data}}(D),
\]

where `alpha_s` and `alpha_W` are the two relative inputs. Before applying
these scales, normalize each precision to a declared mean diagonal or mean
Rayleigh quotient on its effective support. Otherwise equal numerical values
of the two inputs have no comparable meaning.

### 5.4 Why LocScale is not directly another PCG matrix term

LocScale performs windowed amplitude scaling. In each local window, it
multiplies observed Fourier coefficients by a radial scale derived from a
reference-to-observation power ratio while preserving observed phases [1,2].
The scale depends on the current map amplitudes. A direct penalty of the form

\[
\sum_{j,k}\omega_{jk}
\left(\left|F R_j Mx\right|_k-a_{jk}\right)^2
\]

is nonlinear in `x`. Its Hessian is not the fixed symmetric positive-
semidefinite operator required by ordinary PCG. Applying a LocScale transform
inside `apply_normal` would therefore break linearity.

There are two PCG-compatible interpretations:

1. **Wilson covariance route, recommended first.** Use expected local or global
   molecular power to define a fixed covariance/precision. This regularizes
   coefficients but does not attempt literal amplitude replacement.
2. **Lagged LocScale surrogate, later experiment.** At the start of a
   refinement iteration, freeze local complex targets or linearized weights
   from the previous iteration. The resulting quadratic surrogate adds a
   fixed `L_W^T L_W` term and, for a nonzero target, an `L_W^T t_W` right-hand
   side. The targets must remain unchanged throughout both half solves.

The second route can still preserve one solve per half because the surrogate
is prepared before reconstruction. It is nevertheless a nonlinear outer
fixed-point method, costs many overlapping local FFTs, and may anchor phases to
the previous reference. If tested, phase-bearing targets must be generated
independently from the previous even and odd maps. A merged or common
phase-bearing target would directly compromise half-set independence. A common
radial variance curve contains no phase and is the safer first experiment.

### 5.5 Interaction with existing regularization

The Wilson precision partly overlaps with an isotropic Tikhonov or FSC/SSNR
Fourier prior. It must not simply be stacked at an arbitrary strength and then
credited with any improvement. The small base `lambda_0` can remain as a
numerical safeguard, but a Wilson experiment should either disable an optional
FSC/SSNR molecular prior or define the Wilson curve as an explicit replacement
or modulation of that prior.

The solvent precision has less direct overlap because it is spatially confined
to the complement of the molecular envelope. Even so, both priors interact in
the soft transition. Boundary diagnostics and the four-way ablation are
required.

## 6. Refinement workflow

### 6.1 Recommended state/half sequence

For each state:

1. Resolve and validate the previous available state-specific NU envelope.
2. Construct the fixed Wilson profile, covariance approximation, and optional
   mean before either half solve.
3. Reduce or trailing-blend the raw, data-only even statistics.
4. Construct the even PCG operator with the existing broad support, fixed
   solvent weights, and fixed Wilson precision.
5. Finalize and solve the even system once.
6. Repeat the same reduction/finalization/solve sequence once for the odd
   statistics, using the same fixed envelope, Wilson profile, and strengths.
7. Write the normal even and odd maps and their merged map.
8. Calculate FSC/cFAR and resolution metadata from this pair.
9. Run the existing NU filtering and envelope-generation workflow.
10. Make the newly generated envelope available to the next iteration.

This is one solve per half. The solvent feature does not reopen accumulators,
warm-start a second solve, or create a second half-map pair. Existing behavior
of other regularization modes is outside the scope of this proposal and must
not be used as the architectural model for solvent flatness.

### 6.2 Why the NU envelope must lag by one iteration

The current NU policy already supports a lag-one envelope for matching
references. Reusing that artifact avoids an additional ordering cycle. It also
keeps the operator linear: a same-solve envelope update would make the normal
operator depend on the current PCG iterate.

The first iteration therefore remains unchanged when no prior NU envelope is
available. The overall refinement is still a lag-one fixed-point process: the
priored maps can influence the NU envelope generated for the next iteration.
This feedback cannot be removed without an independent reconstruction or a
frozen external envelope. The first experiment must record envelope occupancy
and overlap across iterations so that support shrinkage is visible.

### 6.3 Half-set and bias discipline

Both halves may use the same fixed spatial prior, but that common prior can
correlate their errors. There is no unregularized pair inside this workflow.
Therefore:

- the reported FSC is the FSC of the solvent-regularized estimator;
- an FSC improvement alone is not evidence of additional recovered signal;
- the raw accumulators must remain half-specific and unchanged;
- the envelope must be fixed before either half solve and must not use the
  current half maps;
- common Wilson covariance/profile parameters must contain no phases;
- any phase-bearing molecular target must be lagged and half-specific;
- scientific comparisons must use separate, otherwise identical
  combinations of zero and positive solvent/Wilson strengths;
- comparisons should include map-to-truth or independent map-to-model measures
  when available.

This is the same reason learned or model-derived priors require special care in
half-map validation [10,11].

## 7. Mask contract

The first implementation should accept only the state-specific NU evidence
envelope. It must not use the spherical NU objective support or silently fall
back to `automask3D_stateNN.mrc`.

Before a solve, the owner must validate:

- file presence and complete read;
- box dimensions and sampling distance;
- finite voxel values;
- a documented molecular-envelope convention;
- clipping or rejection outside `[0,1]`;
- nonzero `S` and an adequate effective solvent population;
- compatibility with the current cropped reconstruction grid;
- one fixed mask identity for both halves and all iterations of each solve.

Recommended diagnostics are:

```text
pcg_solvent_prior_enabled=
pcg_solvent_mask_file=
pcg_solvent_mask_min=
pcg_solvent_mask_max=
pcg_solvent_weight_sum=
pcg_solvent_weight_fraction=
pcg_solvent_lambda_rel=
pcg_solvent_lambda_eff=
pcg_solvent_mean_final=
pcg_solvent_rms_final=
pcg_solvent_penalty_final=
pcg_solvent_skip_reason=
pcg_wilson_prior_enabled=
pcg_wilson_mode=
pcg_wilson_lambda_rel=
pcg_wilson_lambda_eff=
pcg_wilson_profile_source=
pcg_wilson_profile_min=
pcg_wilson_profile_max=
pcg_wilson_penalty_final=
```

An invalid or evidence-null mask must produce a clear skip reason. A positive
user setting must never select a different mask type without notice.

## 8. Proposed ownership and change map

| Concern | Owner | Proposed responsibility |
| --- | --- | --- |
| CLI/default/validation | `parameters`, command dictionaries, UI | Define independent solvent and Wilson relative strengths; both default to `0.0`; reject negative values. |
| Numerical priors | `src/main/volume/simple_reconstructor_pcg.f90` | Store fixed solvent and molecular-prior state; add their operator actions in the shared normal-operator wrapper; add a nonzero molecular mean to a separate solve-time RHS. |
| PCG refinement assembly | `src/main/strategies/parallelization/simple_refine3D_strategy.f90` and `simple_rec3D_pcg_strategy.f90` | Resolve one lagged envelope and one fixed Wilson model per state and attach them to the normal even and odd solves without adding a solve. |
| NU envelope generation | `src/main/commanders/simple/simple_commanders_rec_distr.f90` and NU domain modules | Preserve the post-reconstruction generation order and publish the envelope for a later iteration. |
| Workflow policy | `doc/policies/reconstruct3D_pcg_policy.md` and `doc/policies/nonuniform_filtering_policy.md` | Record the single-solve, lag-one-mask, FSC interpretation, missing-mask, and default-off contracts. |

Both precision actions belong behind the high-level `apply_normal` contract,
after the selected data operator is applied and before the outer support
multiplication. This placement gives kernel and matrix-free modes the same
priors without duplicating code or hiding disagreements between the two data
operators.

The raw `(B,D)` file format must not change. The mask and prior strength are
solve-time state and are not accumulated into `D` or persisted in a trailing
chain. A molecular prior mean contributes to a separate solve-time RHS, not to
raw `B`. Loading either prior must not trigger another raw-accumulator read.

## 9. Preconditioning and cost

The exact solvent matvec requires:

- one deterministic double-precision reduction for `s^T x`;
- one multiplication/subtraction pass for `Q_s x`;
- storage for one real solvent-weight volume.

This cost should be small compared with a PCG data-operator FFT. The first
implementation should leave the preconditioner unchanged. This isolates the
effect of the prior and avoids an unjustified approximation of the rank-one
centering term.

A global shell-diagonal Wilson precision needs a forward and inverse volume FFT
unless its action can reuse an existing kernel-operator transform. A local
windowed Wilson or LocScale surrogate needs many small FFTs per matrix-vector
product and is not suitable for the first combined implementation. If a
shell-diagonal Wilson term proves useful, its diagonal spectrum is also a
natural measured candidate for inclusion in the existing FFT preconditioner.

If iteration counts or curvature checks regress, the next step is to assess a
symmetric preconditioner that includes an approximation to `lambda_s D_s`.
The rank-one term must remain in the exact operator. Any preconditioner change
requires its own parity and convergence study.

## 10. Testable hypotheses

- **H1:** At a useful positive strength, the solvent prior reduces
  weighted solvent RMS without a material loss of molecular signal.
- **H2:** With strength zero, outputs are identical to the current route. With
  positive strength, raw `(B,D)` remains unchanged and each half is still
  solved exactly once.
- **H3:** A soft NU envelope is more robust to boundary error than a hard
  solvent projection, especially for weak peripheral density.
- **H4:** The added prior actions have acceptable cost compared with the PCG
  data operator and do not materially increase the iteration count at the
  selected strengths.
- **H5:** Wilson and solvent precisions are complementary: the combined result
  improves molecular spectral behavior and solvent flatness beyond the
  corresponding single-prior runs without additional weak-density loss.

H1, H3, and H5 are scientific hypotheses. H2 is an architectural invariant.
H4 is a performance gate. The half-map FSC is an outcome, not an invariant. Before a
real-data sweep, the experiment owner should record the numerical thresholds
for “material loss,” parity, convergence, and acceptable overhead. This avoids
choosing thresholds after inspecting the maps.

## 11. Staged validation plan

### Gate A: algebra and operator identity

Use small synthetic volumes and fixed masks to verify:

- linearity of `Q_s`;
- the dot-product identity `<x,Q_s y> = <Q_s x,y>`;
- nonnegative quadratic form within floating-point tolerance;
- `Q_s 1 = 0`;
- zero action in voxels with zero solvent confidence;
- continuous changes across a graded mask edge;
- agreement of the implemented penalty gradient with a finite-difference
  derivative of `R_s`;
- unchanged kernel-versus-matrix-free data-operator comparison when the prior
  is tested separately;
- default-off numerical identity with the current PCG route;
- self-adjointness and nonnegative curvature of the Wilson action;
- self-adjointness and nonnegative curvature of the summed prior action;
- correct prior-RHS construction for a synthetic nonzero Wilson mean.

Failure of symmetry, positive curvature, or default-off identity blocks all
workflow tests.

### Gate B: workflow and artifact invariants

Verify that enabling the prior:

- does not change raw `(B,D)` checksums;
- performs exactly one normal solve for each populated half;
- uses the same envelope for the even and odd solves;
- gives shared/distributed parity under the existing deterministic reduction
  tolerance;
- preserves restart and trailing-reconstruction equivalence;
- records a clear skip when the prior envelope is absent;
- generates the current NU envelope only after both solves and FSC;
- reproduces current outputs when both strengths are zero;
- does not add a second solve when either or both priors are enabled.

### Gate C: controlled scientific tests

Start with a synthetic reconstruction with known support and ground truth.
Measure:

- solvent weighted RMS relative to the matched strength-zero run;
- FSC to ground truth, not only half-map FSC;
- molecular-region RMS change;
- boundary ringing and density leakage;
- recovery of a deliberately weak peripheral domain;
- PCG residual history and stop reason.

Repeat with the true envelope, one eroded envelope, one dilated envelope, and
one envelope that omits the weak domain. The omitted-domain case is the main
bias test.

For the first real-data A/B test, sweep

```text
pcg_solvent_lambda_rel =
0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1
```

The range is exploratory, not a proposed production default. Run each point as
a separate one-solve-per-half reconstruction with identical particles, poses,
raw accumulators, PCG iteration limits, and incoming NU envelope. Inspect
solvent variance, weak density, map-model FSC against a held-out or
independently refined model, and convergence. Continue selected strengths for
several refinement iterations and monitor envelope occupancy and overlap to
detect lagged feedback.

After selecting non-destructive individual strengths, run the factorial
ablation:

| Run | Solvent prior | Wilson prior | Purpose |
| --- | --- | --- | --- |
| A | off | off | Current estimator. |
| B | on | off | Isolate solvent flatness. |
| C | off | on | Isolate molecular spectral regularization. |
| D | on | on | Measure complementarity or antagonism. |

All four runs must use the same particles, poses, incoming envelope, Wilson
profile, reconstruction band, and iteration count. Tune neither strength from
half-map FSC alone. Prefer predictive likelihood on held-out particles,
map-to-truth FSC for synthetic data, or an independent map-to-model comparison.
The combined prior passes only if it improves the molecular and solvent
metrics jointly rather than trading one for the other.

### Gate D: acceptance for an experimental release

The feature can remain available as an experimental opt-in only if:

- all algebraic and workflow invariants pass;
- no default-off output changes;
- any FSC increase is supported by an independent truth/model comparison and
  is not treated as sufficient evidence by itself;
- weak-domain recovery is not systematically worse under reasonable mask
  perturbations;
- the combined prior outperforms or lies on a useful Pareto frontier relative
  to both single-prior results;
- the chosen strength reduces solvent variance without material molecular
  signal loss;
- per-half convergence stays within the current operational PCG budget or a
  measured new budget is approved;
- shared and distributed routes agree within the declared tolerance.

No positive production default should be selected from one specimen.

## 12. Main risks and mitigations

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| Envelope misses weak molecular density | Real signal is flattened as solvent. | Use a soft NU envelope, lag it by one iteration, test eroded/omitted masks, and keep the default off. |
| The current priored maps influence the next envelope | Lagged, self-reinforcing support shrinkage. | Track envelope occupancy/overlap, test omitted-domain masks, and add conservative hysteresis only if the first experiment shows shrinkage. |
| Common prior inflates half-map correlation | Misleading resolution estimate. | Compare with a separate strength-zero run and require independent truth/model evidence; do not accept FSC gain alone. |
| Wilson precision duplicates FSC/SSNR or base regularization | Excess molecular shrinkage is attributed to a new physical prior. | Disable the optional overlapping prior in the first ablation or define an explicit replacement/modulation rule. |
| A LocScale amplitude target is treated as linear | The operator changes with `x`, so PCG assumptions fail. | Use a fixed Wilson covariance or a lagged quadratic surrogate; never recompute scaling inside `apply_normal`. |
| A common LocScale target carries phases | Even/odd errors become directly correlated. | Share only phase-free physical profiles; keep any lagged phase-bearing target half-specific. |
| Molecular and solvent precisions have incomparable normalization | The two relative strengths cannot be interpreted. | Normalize each precision on its effective support before multiplying by the common data scale. |
| Outside-support zeros dominate the mean | Prior becomes an unintended zero-density penalty. | Include broad support in `w_v` and exclude negligible support weights. |
| `L` is added instead of `L^T L` | Loss of symmetry and possible PCG failure. | Implement and test `Q_s = D_s - ss^T/S` directly. |
| Solvent population is too small | Unstable mean and ineffective prior. | Require a minimum effective weight and skip with a diagnostic. |
| Large strength worsens conditioning | More iterations or curvature failure. | Start with a log sweep, monitor residuals, and defer preconditioner changes. |
| Silent mask fallback changes the experiment | Results cannot be interpreted. | Permit no fallback in the first implementation. |

## 13. Recommended implementation staging

1. Implement and validate solvent `Q_s` with both strengths zero by default.
2. Add a zero-mean, shell-diagonal Wilson covariance precision with its own
   default-zero strength.
3. Run the four-way ablation with fixed poses and incoming envelope.
4. If the combination is useful, assess whether the Wilson spectrum can enter
   the FFT preconditioner without changing the exact operator.
5. Only then test Singer/Gilles off-diagonal covariance structure.
6. Treat a local LocScale-style quadratic surrogate as a separate research
   feature with half-specific lagged targets and a measured local-FFT cost.

The Wilson covariance and solvent precision are the scientifically clean
combination for the first joint experiment. Literal LocScale amplitude
matching belongs outside the initial linear PCG implementation.

## 14. Recommended reading order

1. Wilson 1942 for the original high-frequency intensity argument [3].
2. Wilson 1949 for the intensity probability distributions [4].
3. Singer 2021 for the modern derivation and cryo-EM interpretation [5].
4. Gilles and Singer 2022 for a Wilson-statistics Bayesian prior [6].
5. Wang 1985 for the original practical solvent-flattening method [7].
6. Terwilliger 1999 for the explicit likelihood view used by this proposal
   [9].
7. Terwilliger et al. 2020 for the cryo-EM density-modification context [10].
8. LocScale-2.0 for the immediate motivation and its stated real-space-prior
   direction [1].

## References

1. A. Bharadwaj, R. de Bruin, and A. J. Jakobi, “Confidence-guided cryo-EM
   map optimisation with LocScale-2.0,” *Nature Communications*, vol. 17,
   article 8778, 2026. [doi:10.1038/s41467-026-75327-8](https://doi.org/10.1038/s41467-026-75327-8)
   ([open-access article](https://www.nature.com/articles/s41467-026-75327-8)).

2. A. J. Jakobi, M. Wilmanns, and C. Sachse, “Model-based local density
   sharpening of cryo-EM maps,” *eLife*, vol. 6, e27131, 2017.
   [doi:10.7554/eLife.27131](https://doi.org/10.7554/eLife.27131)
   ([open-access article](https://elifesciences.org/articles/27131)).

3. A. J. C. Wilson, “Determination of absolute from relative X-ray intensity
   data,” *Nature*, vol. 150, p. 152, 1942.
   [doi:10.1038/150152a0](https://doi.org/10.1038/150152a0).

4. A. J. C. Wilson, “The probability distribution of X-ray intensities,”
   *Acta Crystallographica*, vol. 2, pp. 318-321, 1949.
   [doi:10.1107/S0365110X49000813](https://doi.org/10.1107/S0365110X49000813).

5. A. Singer, “Wilson statistics: derivation, generalization and applications
   to electron cryomicroscopy,” *Acta Crystallographica Section A*, vol. 77,
   pp. 472-479, 2021.
   [doi:10.1107/S205327332100752X](https://doi.org/10.1107/S205327332100752X)
   ([open-access full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8477642/)).

6. M. A. Gilles and A. Singer, “A molecular prior distribution for Bayesian
   inference based on Wilson statistics,” *Computer Methods and Programs in
   Biomedicine*, vol. 221, article 106830, 2022.
   [doi:10.1016/j.cmpb.2022.106830](https://doi.org/10.1016/j.cmpb.2022.106830)
   ([open-access full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9233040/)).

7. B.-C. Wang, “Resolution of phase ambiguity in macromolecular
   crystallography,” *Methods in Enzymology*, vol. 115, pp. 90-112, 1985.
   [doi:10.1016/0076-6879(85)15009-3](https://doi.org/10.1016/0076-6879(85)15009-3).

8. K. D. Cowtan and P. Main, “Improvement of macromolecular electron-density
   maps by the simultaneous application of real and reciprocal space
   constraints,” *Acta Crystallographica Section D*, vol. 49, pp. 148-157,
   1993.
   [doi:10.1107/S0907444992007698](https://doi.org/10.1107/S0907444992007698).

9. T. C. Terwilliger, “Reciprocal-space solvent flattening,” *Acta
   Crystallographica Section D*, vol. 55, pp. 1863-1871, 1999.
   [doi:10.1107/S0907444999010033](https://doi.org/10.1107/S0907444999010033)
   ([open-access article](https://journals.iucr.org/d/issues/1999/11/00/gr0940/)).

10. T. C. Terwilliger, O. V. Sobolev, P. V. Afonine, P. D. Adams, and
    R. J. Read, “Density modification of cryo-EM maps,” *Acta
    Crystallographica Section D*, vol. 76, pp. 912-925, 2020.
    [doi:10.1107/S205979832001061X](https://doi.org/10.1107/S205979832001061X)
    ([open-access article](https://journals.iucr.org/d/issues/2020/10/00/ir5011/)).

11. D. Kimanius, G. Zickert, T. Nakane, J. Adler, S. Lunz,
    C.-B. Schönlieb, O. Öktem, and S. H. W. Scheres, “Exploiting prior
    knowledge about biological macromolecules in cryo-EM structure
    determination,” *IUCrJ*, vol. 8, pp. 60-75, 2021.
    [doi:10.1107/S2052252520014384](https://doi.org/10.1107/S2052252520014384)
    ([open-access article](https://journals.iucr.org/m/issues/2021/01/00/fq5015/)).
