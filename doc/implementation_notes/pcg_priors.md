# PCG priors: experiment record and the solvent/Wilson proposal

## Status

This is the single reference for regularization of the PCG reconstruction
backend: what exists on master, what was tried and removed (with the
measurements that killed it), and the approved-for-staging proposal for a
solvent-flatness prior and a Wilson molecular prior. It consolidates and
replaces `pcg_euclid_crash_investigation.md`,
`pcg_reconstruction_production_readiness.md`,
`pcg_refine3D_integration_plan.md`, and
`pcg_real_space_solvent_flatness_prior_proposal.md`. Current implemented
contracts live in `doc/policies/reconstruct3D_pcg_policy.md` and the PCG
section of `doc/policies/refine3D_policy.md`.

**All proposed features are opt-in and disabled by default.** Neither
hard-zeros the solvent region.

## 1. Current regularization inventory (master)

The production PCG path performs, per state and half, one *base* solve and one
*ML replay* solve from the same raw `(B,D)` statistics:

- **base solve**: `(H_data + lambda_0 I) x = b`, cold-started; its half pair
  (`*_unfil`) produces the FSC.
- **ML replay**: `(H_data + P_tau + lambda_0 I) x = b`, where `P_tau` is the
  FSC/SSNR shell-diagonal Fourier precision (`ml_prior`, present in both the
  operator and the preconditioner). Its output is the shipped map. It
  warm-starts from the previous refinement iteration's ML half map when one
  exists (own half only, support re-masked after constant-FOV resampling,
  first-iteration `startvol` excluded), falling back to the base solution.
- `lambda_0` is the fixed absolute base Tikhonov `PCG_LAMBDA = 1e-3`.
- the soft spherical support `P` constrains the solve as `P H P u = P b`.
- maps are stored at the data-quotient amplitude convention
  (`drop_legacy_box_division.md`); deapodization is inside the solver.
- `report_beyond_band_excess` logs when post-band RMS exceeds 10x the
  band-edge shell; it is the regression signal for beyond-band defects.

## 2. Record of prior experiments (tried and removed)

Kept so no experiment is re-run. Measurements in
`drop_legacy_box_division.md` §11–12.

| Experiment | Outcome | Lesson |
| --- | --- | --- |
| `pcg_lambda_lap` — biharmonic (nu^4) smoothness prior, operator + preconditioner | Removed. On the synthetic fixture it suppressed beyond-band noise 10–80x and closed the 2-iteration in-band deficit (ratio 0.84 → 1.00). On real data (bgal, box 256) no usable strength exists: rho falls ~3 orders across the band while nu^4 spans ~1.5 from half-band to edge, so any lambda that suppresses beyond-band also attenuates the band edge (−25% at 5.4 A even at 0.01). | A data-scale-anchored spectral prior cannot separate in-band from beyond-band; band-selective suppression is the FSC-informed `P_tau`'s job. The real benefit was *conditioning* the ML replay (RESID 6.5 → 0.86 at 2 its) — captured instead by the cross-iteration warm start. |
| `pcg_lambda_rel` — relative base-lambda CLI (the "P0 lambda contract") | CLI removed as unused; the internal `set_lambda_relative` mechanism and the deterministic `s_data(D)` data-scale functional remain (test commanders, diagnostics). | Relative anchoring machinery is retained where priors need it; it is not user-facing. |
| Kernel-derived preconditioner | Removed. Did not improve runtime at one iteration and degraded reconstruction metrics. | The sampling-density diagonal with a shell-relative floor is the preconditioner. |
| Whole-volume-L2 ML-vs-base rescale guard | Removed. Whole-volume norms are solvent-dominated and never fired on the failing solves. | Compare shell profiles, not whole-map statistics. |
| `x0 = P b` base-solve warm start (proposed, never built) | Superseded. | The cross-iteration ML warm start addresses the same under-convergence where it matters (the shipped ML map). |
| `ref_taper` — legacy KB matching-reference taper | Removed after a master reliability control (1/10 failures on master vs 2/10 on the branch) showed no regression to fix. | Measure the baseline before attributing a failure rate to a change. |

Two historical defects, both fixed, whose lessons are contracts now:

- **Beyond-band excess (defect A).** PCG base solves can leave hot spectral
  content just above the matching band; a stage transition that extends `kto`
  then exposes it to euclid matching. The diagnostic survives
  (`PCG BEYOND-BAND EXCESS`); the structural mitigation is ML-replay
  convergence (warm start + adequate iterations), not a smoothness prior.
- **Amplitude-convention discontinuity (defect B).** The euclid/sigma2
  equilibrium survives any *stable* reference amplitude but not an abrupt
  scale change between consecutive iterations. Retired by the shared
  data-quotient convention; any future backend or prior must preserve
  scale continuity at handoffs.

Prior priority (updated from the retired production-readiness note):

| Idea | Priority | Recommendation |
| --- | --- | --- |
| solvent-flatness quadratic prior | first new-prior experiment | this proposal, §4 |
| Wilson molecular precision | with solvent, as a `P_tau` spectrum source | this proposal, §5 |
| soft state weights | integration semantics | weight both `B` and `D` |
| symmetry projection/permutation | performance only | after distributed profiling |
| total variation / non-negativity / wavelet L1 | deferred | separate nonlinear/proximal solver project |
| joint pose/shift refinement | deferred | separate program (`continuous_3D_refinement_on_pcg_operator.md`) |

Any prior must preserve half independence, positive-semidefinite CG structure,
conditioning, data-mass scaling, and separate reporting of fit and prior
energy, and must be evaluated against a strength-zero control.

## 3. Decision requested

Approve a staged experiment with the following fixed design:

1. Add a dimensionless control named `pcg_solvent_lambda_rel` with default
   value `0.0`.
2. Interpret a positive value as the strength of a graded solvent-flatness
   prior relative to the PCG data scale.
3. Use the previous available `nu_envmask3D_stateNN.mrc` as a molecular
   envelope. Use its complement as solvent confidence.
4. Hold the envelope fixed during each PCG solve and use the same envelope for
   both half sets.
5. **Add the solvent precision to the ML replay solves only.** The base
   solves — and therefore the FSC — remain unregularized by the new priors.
   Do not add another solve or raw-accumulator reload for this feature; the
   existing cross-iteration ML warm start is retained unchanged.
6. Calculate the ordinary FSC from the base (`_unfil`) pair as today, then run
   the existing NU filtering and envelope-generation steps on the priored ML
   pair.
7. Treat the envelope generated in the current iteration as input to a later
   iteration, never as mutable state inside the current PCG solve.
8. If the NU envelope has a different grid but the same physical extent,
   resample it by the factor-free constant-FOV Fourier pad/clip and re-clip
   values to `[0,1]`; skip the priors for that iteration only when the
   envelope is missing, invalid, or of mismatched physical extent. Do not
   silently substitute a density mask or a spherical mask.
9. Implement the Wilson molecular precision as a second *spectrum source* for
   the existing `P_tau` mechanism (see §5), controlled by
   `pcg_wilson_lambda_rel`, default `0.0` — not as a parallel prior array.
10. Test four separate cases: neither prior, solvent only, Wilson only, and
    both. Do not infer the value of the combination from either prior alone.

Two verification prerequisites before any implementation:

- **The `P H P` contract.** The algebra below assumes the support operator is
  applied on both sides of the normal operator and once on the RHS
  (`reconstruct3D_pcg_policy.md` §4 documents this as the contract); verify
  the shared `apply_normal` wrapper and RHS construction actually satisfy it
  before building on it.
- **The deapodization domain.** With `l_deapod` off, the kernel path wraps the
  KB envelope around the operator and the iterate lives in the enveloped
  domain, where real-space solvent weights defined on the map domain would be
  wrong. Restrict the feature to the deapodized kernel mode (hard error
  otherwise), as `pcg_lambda_lap` did.

## 4. The solvent-flatness prior

### 4.1 Motivation and literature position

LocScale and solvent flattening solve related but different problems.
LocScale-1.0 performs local radial Fourier-amplitude scaling against a model
map [2]. The physics-informed route in LocScale-2.0 constructs a pseudoatomic
or hybrid reference and scales local Fourier amplitudes while preserving the
observed local phases [1]; it lists solvent flattening as a possible *future*
real-space prior. Wilson statistics [3–6] underlies the molecular prior of §5.
Solvent flattening proper starts with Wang [7], is combined with
reciprocal-space constraints by Cowtan and Main [8], expressed as an explicit
likelihood by Terwilliger [9], and carried into cryo-EM density modification
[10]. That likelihood line is the direct ancestor of this penalty. The
proposal keeps three ideas separate:

- **Wilson prior:** a molecular mean/covariance model.
- **LocScale operation:** local, phase-preserving amplitude rescaling against
  an expected or reference spectrum (not PCG-compatible as-is, §5.4).
- **Solvent prior:** expected low real-space variance in the solvent.

### 4.2 Graded solvent confidence

Let `m_v` be the NU molecular envelope at voxel `v` (`0 <= m_v <= 1`,
`1` = confident molecular region), and `p_v` the existing broad support value.
Define

```text
w_v = p_v (1 - m_v),   s_v = w_v^2,   S = sum_v s_v
```

Including `p_v` prevents the many fixed zero voxels outside the broad support
from dominating the solvent mean, and grades the prior smoothly at both the
NU-envelope boundary and the spherical-support boundary.

### 4.3 Weighted solvent mean and penalty

For a trial map `x`, the weighted solvent mean and penalty are

```text
mu_s(x) = (sum_v s_v x_v) / S
R_s(x)  = (lambda_s/2) sum_v s_v [x_v - mu_s(x)]^2
```

The generating operator is `L = W C_s` with `W = diag(w)` and
`C_s = I - 1 s^T / S` (weighted centering). **The normal operator must add
`L^T L`, not `L`**: `L` is not symmetric because weighted centering is not an
orthogonal projection in the ordinary inner product, and adding it directly
would invalidate CG's symmetry assumption. The required precision is

```text
Q_s = L^T L = D_s - s s^T / S,   D_s = diag(s)
Q_s x = s .* [x - mu_s(x)]
```

No matrix is stored: one deterministic weighted reduction and two volume
passes per application.

### 4.4 Per-half ML-replay normal equation

The ML replay for each half becomes

```text
P (H_data + P_tau + lambda_0 I + lambda_s Q_s) P u = P b
lambda_s = pcg_solvent_lambda_rel * s_data(D)
```

The base solves are unchanged. The data-scale anchoring uses the retained
internal `s_data(D)` functional. Note this anchoring is defensible here where
it was not for the removed nu^4 prior (§2): the solvent penalty competes with
the data term as a real-space aggregate over a fixed region, not per Fourier
shell, so the orders-of-magnitude shell-wise decline of rho does not create an
empty strength window. Before applying the scale, normalize the precision to a
declared mean diagonal on its effective support so the two relative inputs of
§5 have comparable meaning.

### 4.5 Algebraic properties

For every real `x`: `x^T Q_s x = sum_v s_v [x_v - mu_s(x)]^2 >= 0`, so `Q_s`
is symmetric positive semidefinite, and `Q_s (c 1) = 0` — it penalizes solvent
*variation*, not the solvent offset, which cryo-EM data do not determine. The
broad support `P` remains: the solvent operator has an unpenalized constant
mode and the interpolation/deapodize path is least reliable near the padded
corners, so removing `P` would add weakly constrained unknowns while the new
prior is being evaluated.

A free hypothesis this null space suggests (H6, §7): the known PCG low-shell
anomalies — the k=1–2 backend amplitude excess (~1.3–2.0) and the centre-bin
radial deficit — are plausibly solvent-offset/low-k noise behavior that `Q_s`
addresses directly. The gated `rec3D_backends` shell and centre-bin
diagnostics measure this for free. Re-baseline those numbers after the ML
warm start is validated, since it moves them too.

## 5. The Wilson molecular prior

### 5.1 Combined model

The combination of solvent and molecular priors is a product-of-experts: the
likelihood and both priors contribute additive negative log probabilities in
complementary regions, and both remain fixed during a half solve, so the
summed precision stays symmetric positive semidefinite. With
`M = diag(m)` the soft molecular selector and `Q_W` a fixed Wilson-derived
precision, the combined ML-replay system is

```text
P [H_data + P_tau + lambda_0 I + lambda_s Q_s + lambda_W M Q_W M] P u
    = P [b + lambda_W M Q_W mu_W]
```

For a zero Wilson mean (`mu_W = 0`, the recommended first form) the RHS is the
data-only `b`. The envelope divides responsibility softly: `M` localizes the
molecular prior, `w = p(1-m)` localizes solvent flatness; both act with
reduced weight in the transition zone, deliberately, because a hard partition
would ring at one mask contour.

### 5.2 Implement Wilson as a `P_tau` spectrum source

A shell-diagonal Wilson precision
`Q_W = F* diag(q_W(k)) F`, with `q_W(k)` the regularized inverse expected
molecular variance, is *structurally the same object* as the existing
FSC/SSNR `ml_prior` — a per-shell Fourier diagonal in both the operator and
the preconditioner. Do not build a second array with its own plumbing:
implement Wilson as an alternative or modulating **spectrum source** for the
one existing `P_tau` mechanism. This settles §5.5's overlap problem by
construction (there is one shell-diagonal precision, with a declared spectrum
provenance: FSC-derived, Wilson-derived, or a documented combination), puts
Wilson into the preconditioner for free, and follows the repository's
one-path rule.

Singer derives the non-diagonal covariance of Fourier coefficients from the
random bag-of-atoms model [5]; Gilles and Singer use the corresponding mean
and covariance as a Bayesian molecular prior [6]. The shell-diagonal
approximation is the lowest-risk first form; off-diagonal structure [5,6] and
LocScale-informed *local* covariance profiles (sums of
`M R_j^T F* diag(q_W,j) F R_j M`, each positive semidefinite) are later
experiments — the local form costs many small FFTs per matvec.

A zero-mean choice adds only precision — it suppresses coefficients according
to expected variance and does **not** force a target amplitude spectrum. A
nonzero Wilson mean adds a fixed prior term to the RHS and needs stronger
provenance and bias controls.

### 5.3 Why LocScale is not directly another PCG matrix term

LocScale multiplies observed local Fourier amplitudes by a scale derived from
a reference-to-observation power ratio while preserving phases [1,2]; the
scale depends on the current map, so the corresponding penalty is nonlinear in
`x` and its Hessian is not a fixed SPD operator — applying it inside
`apply_normal` breaks linearity. The PCG-compatible interpretations are (a)
the fixed Wilson covariance above (recommended first) and (b) a *lagged*
quadratic surrogate frozen at iteration start — a nonlinear outer fixed-point
method with many local FFTs, whose phase-bearing targets must be half-specific
and lagged (a merged phase-bearing target would directly compromise half-set
independence; a common radial variance curve contains no phase and is the
safer first experiment).

## 6. Workflow, mask contract, and bias discipline

### 6.1 Per-state sequence (one particle accumulation, two solve phases)

1. Resolve and validate (or constant-FOV-resample, §3 item 8) the previous
   state-specific NU envelope; construct the fixed Wilson spectrum before
   either half solve.
2. Reduce the raw, data-only even and odd statistics as today.
3. Base-solve each half unchanged (`P_tau = 0`, no new priors), write the
   `_unfil` pair, and compute FSC/cFAR from it — the FSC keeps its current,
   unregularized meaning.
4. Build `P_tau` from the FSC (or the Wilson spectrum source, per §5.2), add
   `lambda_s Q_s`, and run each half's ML replay with the existing
   cross-iteration warm start.
5. Write the standard maps, run NU filtering and envelope generation, and
   publish the new envelope for a later iteration.

### 6.2 Mask contract

Accept only the state-specific NU evidence envelope; never fall back to the
spherical NU support or `automask3D_stateNN.mrc` silently. Validate before a
solve: file presence and complete read; finite values; documented
molecular-envelope convention; clip/reject outside `[0,1]`; nonzero `S` and an
adequate effective solvent population; same physical extent as the
reconstruction grid (resample per §3 item 8 when the lattice differs); one
fixed mask identity for both halves and all iterations of each solve. Emit the
diagnostics:

```text
pcg_solvent_prior_enabled=  pcg_solvent_mask_file=   pcg_solvent_mask_min/max=
pcg_solvent_weight_sum=     pcg_solvent_weight_fraction=
pcg_solvent_lambda_rel/eff= pcg_solvent_mean_final=  pcg_solvent_rms_final=
pcg_solvent_penalty_final=  pcg_solvent_skip_reason=
pcg_wilson_prior_enabled=   pcg_wilson_mode=         pcg_wilson_lambda_rel/eff=
pcg_wilson_profile_source=  pcg_wilson_profile_min/max= pcg_wilson_penalty_final=
```

An invalid or evidence-null mask must produce a clear skip reason; a positive
user setting must never select a different mask type without notice.

### 6.3 Half-set and bias discipline

Both halves may use the same fixed spatial prior, but a common prior can
correlate their errors. With the priors confined to the ML replay, the
reported FSC remains that of the unregularized base pair — but the lag-one
feedback through the NU envelope remains: priored maps influence the next
iteration's envelope. Therefore: record envelope occupancy/overlap across
iterations so support shrinkage is visible; keep raw accumulators
half-specific and unchanged; fix the envelope before either half solve; share
only phase-free Wilson parameters; make scientific comparisons against
separate strength-zero runs; and prefer map-to-truth or independent
map-to-model measures over half-map FSC [10,11].

### 6.4 Ownership

| Concern | Owner |
| --- | --- |
| CLI/defaults/validation | `parameters`, dictionaries, UI (both strengths default `0.0`, reject negatives) |
| Numerical priors | `src/main/volume/simple_reconstructor_pcg.f90` (fixed solvent weights and Wilson spectrum state; actions behind the shared normal-operator wrapper, after the data operator and inside the support sandwich) |
| PCG refinement assembly | `simple_rec3D_pcg_strategy.f90` (resolve one lagged envelope and one fixed spectrum per state; attach to the ML replays) |
| NU envelope generation | `simple_commanders_rec_distr.f90` and NU modules (unchanged order; publish for a later iteration) |
| Policy record | `reconstruct3D_pcg_policy.md`, `nonuniform_filtering_policy.md` |

The raw `(B,D)` format does not change; mask and strengths are solve-time
state, never accumulated into `D`; a nonzero Wilson mean contributes to a
solve-time RHS, not to raw `B`; loading a prior triggers no raw re-read.

## 7. Hypotheses

- **H1:** At a useful positive strength, the solvent prior reduces weighted
  solvent RMS without material loss of molecular signal.
- **H2 (architectural invariant):** With both strengths zero, outputs are
  identical to the current route; with positive strengths, raw `(B,D)`, the
  base solves, and the FSC are unchanged, and the solve structure remains one
  base solve plus one ML replay per half.
- **H3:** A soft NU envelope is more robust to boundary error than a hard
  solvent projection, especially for weak peripheral density.
- **H4 (performance gate):** Prior actions cost little against the PCG data
  operator, and the production iteration budget still converges the ML replay
  (measured, not assumed — see Gate C).
- **H5:** Wilson and solvent precisions are complementary: the combination
  improves molecular spectral behavior and solvent flatness jointly, not by
  trading one for the other.
- **H6:** The solvent prior reduces the PCG k=1–3 amplitude excess and the
  centre-bin deficit measured by the gated `rec3D_backends` diagnostics.

Record numerical thresholds for "material loss", parity, convergence, and
acceptable overhead *before* the real-data sweep.

## 8. Staged validation

### Gate A: algebra and operator identity

Small synthetic volumes, fixed masks: linearity of `Q_s`; the dot-product
identity `<x, Q_s y> = <Q_s x, y>`; nonnegative quadratic form; `Q_s 1 = 0`;
zero action where solvent confidence is zero; continuity across a graded mask
edge; penalty gradient vs finite differences; unchanged
kernel-vs-matrix-free data-operator comparison; default-off numerical
identity; self-adjointness and nonnegative curvature of the Wilson action and
of the summed prior action; correct prior-RHS for a synthetic nonzero Wilson
mean. Failure of symmetry, curvature, or default-off identity blocks all
workflow tests.

### Gate B: workflow and artifact invariants

Enabling the priors: does not change raw `(B,D)` checksums; leaves the base
solves and the FSC bit-identical; performs exactly one ML replay per populated
half; uses the same envelope for even and odd; keeps shared/distributed
parity under the existing deterministic-reduction tolerance; records a clear
skip when the envelope is absent; generates the NU envelope only after both
solve phases and the FSC; reproduces current outputs at strength zero.
(Fractional/trailing reconstruction is hard-errored on the PCG path today,
which keeps this gate's surface small; revisit when that guard lifts.)

### Gate C: controlled scientific tests — with convergence isolation

The direct lesson of the removed nu^4 prior: at production iteration budgets,
measured "prior effects" are confounded with ML-replay under-convergence, and
`D_s` is real-space diagonal — not representable in the Fourier-diagonal
preconditioner — so some per-iteration convergence cost is guaranteed.
**Scientific comparisons run at converged settings** (high `maxits_pcg` /
tight `rtol`); the production-budget behavior is measured separately as H4.

Harness: the neutral-phantom fixture and the gated `rec3D_backends`
ground-truth mode (map-to-truth FSC, radial LS profiles, background handling,
centre-bin and shell diagnostics) already provide the measurement
infrastructure; the eroded/dilated/omitted-envelope runs drop straight into
it. Measure: weighted solvent RMS vs the matched strength-zero run; FSC to
ground truth (not only half-map FSC); molecular-region RMS change; boundary
ringing and leakage; recovery of a deliberately weak peripheral domain
(the omitted-domain envelope is the main bias test); PCG residual history and
stop reason; the H6 low-shell/centre diagnostics.

First real-data sweep (exploratory, not a production default):

```text
pcg_solvent_lambda_rel = 0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1
```

each point a separate reconstruction with identical particles, poses, raw
accumulators, converged solver settings, and incoming NU envelope. Then the
factorial ablation:

| Run | Solvent | Wilson | Purpose |
| --- | --- | --- | --- |
| A | off | off | Current estimator. |
| B | on | off | Isolate solvent flatness. |
| C | off | on | Isolate molecular spectral regularization. |
| D | on | on | Complementarity or antagonism. |

All four with the same particles, poses, incoming envelope, Wilson spectrum,
band, and solver settings. Tune neither strength from half-map FSC alone;
prefer held-out predictive likelihood, map-to-truth FSC, or independent
map-to-model comparison. Continue selected strengths across refinement
iterations and monitor envelope occupancy/overlap for lagged feedback.

### Gate D: acceptance for an experimental release

All algebraic and workflow invariants pass; no default-off output changes; FSC
increases are supported by independent truth/model evidence; weak-domain
recovery is not systematically worse under reasonable mask perturbations; the
combination outperforms or sits on a useful Pareto frontier vs both single
priors; the chosen strength reduces solvent variance without material
molecular loss; per-half convergence stays within the operational budget or a
measured new budget is approved; shared and distributed routes agree. No
positive production default from one specimen.

## 9. Risks and mitigations

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| Envelope misses weak molecular density | Real signal flattened as solvent. | Soft NU envelope, lag-one, eroded/omitted-mask tests, default off. |
| Priored maps influence the next envelope | Lagged self-reinforcing support shrinkage. | Track occupancy/overlap; omitted-domain tests; hysteresis only if shrinkage appears. |
| Common prior inflates half-map correlation | Misleading resolution estimate. | Priors in the ML replay only (FSC stays base-pair); independent truth/model evidence required. |
| Wilson duplicates the FSC/SSNR precision | Shrinkage credited to a new physical prior. | One `P_tau` mechanism with a declared spectrum source (§5.2); ablation. |
| LocScale amplitude target treated as linear | Operator depends on `x`; PCG assumptions fail. | Fixed Wilson covariance or a lagged surrogate; never rescale inside `apply_normal`. |
| Common LocScale target carries phases | Half errors directly correlated. | Share only phase-free profiles; lagged, half-specific otherwise. |
| Incomparable prior normalizations | Relative strengths uninterpretable. | Normalize each precision on its effective support before the common data scale. |
| Outside-support zeros dominate the mean | Unintended zero-density penalty. | Include broad support in `w_v`. |
| `L` added instead of `L^T L` | Loss of symmetry; PCG failure. | Implement `Q_s = D_s - s s^T/S` directly; Gate A identity. |
| Solvent population too small | Unstable mean, ineffective prior. | Minimum effective weight; skip with diagnostic. |
| Prior effect confounded with under-convergence | Wrong scientific conclusions (the nu^4 lesson). | Gate C convergence isolation; H4 measured separately. |
| Preconditioner cannot represent `D_s` | Iteration-count regression at production budgets. | Measure; if material, assess a symmetric approximation to `lambda_s D_s` (rank-one term stays exact) with its own parity study. |
| Envelope grid mismatch at stage handoffs | Prior silently off for a stage. | Constant-FOV resample + `[0,1]` re-clip (§3 item 8). |
| Silent mask fallback | Uninterpretable experiment. | No fallback in the first implementation. |

## 10. Staged development plan

Each stage ends at a gate; a later stage does not begin until the preceding
gate passes and the stage's changes are committed (trunk development;
everything default-off, so master is protected throughout). "Fixture" means
the neutral-phantom project with ground truth plus the gated
`test=rec3D_backends`.

### Cross-cutting rules (the distilled foot-gun list)

- **R1 — one variable per experiment.** Never change solver behavior and a
  baseline in the same measurement.
- **R2 — no moving baselines.** Re-record every reference number after any
  solver change lands (the warm start moves the centre-bin, k=1–3, and RESID
  numbers the priors will be judged by).
- **R3 — converged-settings science, production-budget cost.** Scientific
  claims at high `maxits_pcg`/tight `rtol`; iteration-budget impact measured
  separately (the nu^4 lesson: 2-iteration results inverted conclusions).
- **R4 — every gate gets a mutation test where feasible.** A gate that has
  never failed on a deliberate defect is not known to gate anything
  (the box-factor/deapodization mutations are the model).
- **R5 — strength-zero bit-identity at every stage boundary.** Raw `(B,D)`
  checksums, `_unfil` maps, and FSC unchanged with the feature off.
- **R6 — never tune from half-map FSC**; truth/model comparisons or held-out
  evidence decide.
- **R7 — reliability changes require a measured control** at matched n (the
  `ref_taper` lesson: master's own failure rate was 1/10).
- **R8 — plumbing before math.** Ship inert loaders/diagnostics first so mask
  handling and operator algebra are never debugged simultaneously.
- **R9 — record acceptance thresholds in this document before running the
  experiment that is judged by them.**
- **R10 — every prior is an independent toggle.** One strength control per
  prior, default `0.0` = off, forwarded intact through `abinitio3D` ->
  `refine3D` -> the solver, with any combination selectable. This is both the
  experimental design (factorial ablation by command line) and the
  application path: refine3D is the work-horse everything else depends on,
  so a validated prior is deployed by turning its flag on — no new
  integration step.

### Stage 0 — close out in-flight work (blocking)

The cross-iteration ML warm start is implemented and unvalidated; it changes
the baselines (R2).

- 0.1 **VALIDATED (2026-08-25):** the cross-iteration warm start produces a
  clean bgal map with a sensible resolution estimate at production settings.
- 0.2 Re-baseline and record here: fixture, streptavidin, and bgal
  `rec3D_backends` numbers (in-band ratio, k=1–3, centre-bin, RESID) at both
  2 and converged iterations.
- 0.3 Commit the consolidation + warm start.

### Stage 1 — verification prerequisites (no behavior change)

- 1.1 **DONE (2026-08-25).** Audit confirmed: `apply_normal` applies the
  support on both sides and the RHS carries exactly one `P` in both
  `end_accum` and `solve`. `test=pcg_recon` Stage 3b now asserts the masked
  operator's dot-product symmetry and strict positivity (the soft edge makes
  `P` non-idempotent, so a one-sided mask application fails this gate where
  the unmasked check cannot see it). Full suite green.
- 1.2 **DONE (2026-08-25).** Production never calls `set_deapod` (default on)
  and always sets `PCG_OP_KERNEL`. `assert_prior_attachment_mode` (hard error
  unless kernel + deapodized) is called at both ML-replay attach sites —
  shared `regularize_state_half` and distributed `reduce_solve_state_half`.
- 1.3 **DONE, refine3D diagnostics check pending (2026-08-25).** Inert
  `report_solvent_envelope_status` in the strategy runs the S6.2 contract
  (presence, cubic lattice, physical-extent identity, constant-FOV
  `read_and_crop` resample on lattice mismatch, finiteness, `[0,1]` re-clip,
  nonzero solvent evidence) and emits the `pcg_solvent_*` block with
  `prior_enabled=F`, once per state before either half's ML replay, both
  execution paths. Fixture gate PASS unchanged (median ratio .8428, identical
  to the pre-change runs). NOTE: the fixture runs `objfun=cc`, and
  `l_ml_reg` requires `objfun=euclid`, so the diagnostic block cannot fire
  there — the envelope-absent skip and the resample demo must be read from
  the next production refine3D/abinitio3D run.

### Stage 2 — solvent `Q_s` operator, Gate A (algebra), default off

**DONE (2026-08-25).** `Q_s` lives in `simple_reconstructor_pcg`
(`set_solvent_prior` / `apply_solvent_precision` / `get_solvent_stats`):
matrix-free `s.*(x - mu_s(x))` with `s = [p(1-m)]^2` normalized to unit mean
diagonal on its effective support; effective strength
`lambda_solvent = pcg_solvent_lambda_rel * s_data(D)` derived alongside the
relative ridge lambda in `update_lambda_from_density`; attached in both
concrete operators (deliberately NOT in the Fourier-shell preconditioner — a
real-space diagonal cannot be fused there; conditioning, not correctness).
`pcg_solvent_lambda_rel` (default 0.0) is registered and exposed on
`reconstruct3D`; the strategy attaches the prior to both halves' ML replays
when the envelope loader (`resolve_solvent_envelope`, the Stage-1.3 loader
gone live) validates it AND the strength is positive, with a per-half
`PCG SOLVENT PRIOR` line (`lambda_eff`, `mean/rms/penalty_final`) after each
priored solve, and a loud THROW_WARN skip when a positive strength has no
usable envelope.

- Gate results: `test=pcg_priors` (8 stages: normalization contract, adjoint,
  PSD, `Q_s 1 = 0`, zero action on the `m=1` plateau, graded-edge continuity
  at eps=1e-3, FD gradient of `R_s`, composition with the masked normal
  operator) ALL PASS; mutation `L` instead of `L^T L` (outer weight
  `sqrt(s)`) fails the adjoint, gradient, and composition stages; mutation
  dropping the rank-one mean term fails exactly the null-space stage (and
  only it — `s.*x` is still symmetric PSD, the correct signature);
  `test=pcg_recon` full suite unchanged; R5 fixture bit-identity holds
  (gated median .8428, identical to baseline, strength defaulting to 0).
- Remaining for Stage 3: shared vs `nparts=2` parity at a positive strength,
  abinitio3D/refine3D forwarding, `rec3D_backends` strength registration.

### Stage 3 — workflow integration incl. abinitio3D, Gate B (invariants)

**Implementation DONE (2026-08-25); real-data gate items pending.**
`pcg_solvent_lambda_rel` is registered on `refine3D` and `abinitio3D`
(activation-gated on `rec_backend=pcg`, like `maxits_pcg`) and on
`test=rec3D_backends`; abinitio3D forwards it through the inherited
`cline_refine3D` plus `apply_refine3D_reconstruction_controls` (R10, exactly
the `maxits_pcg`/`rtol` route). No further forwarding was needed: the PCG
master always runs inside the process whose own command line carries the key
(`simple_refine3D_strategy` and `simple_rec3D_strategy` both pass `params`),
and workers only accumulate raw statistics.

- Gate results so far: `test=pcg_priors` Stage 9 solves the same priored
  synthetic problem through monolithic streaming AND two-part raw reduction
  — the exact seam separating shared from `nparts=2` execution — with
  identical effective strengths and rel_err(x) = 4.2e-5 (PASS); the fixture
  `rec3D_backends` with `pcg_solvent_lambda_rel=0` explicit on the command
  line is bit-identical (.8428) — R5 strength-zero identity end-to-end.
- **Pending on real data** (user-side runs): envelope-absent skip lines in a
  production refine3D log; one full abinitio3D run with
  `rec_backend=pcg automsk=yes` and a positive strength completing, showing
  the prior activating exactly when the envelope first exists lag-one and
  exercising the stage-handoff resample; gated `rec3D_backends` at a small
  positive strength (needs `objfun=euclid` + envelope, i.e. a refine3D run
  directory).

### Stage 4 — Gate C science, synthetic (converged settings, R3)

- Record thresholds first (R9): "material molecular loss", acceptable
  weighted-solvent-RMS reduction, omitted-domain damage bound.
- Fixture sweep at converged settings; envelope variants: true, eroded,
  dilated, omitted-domain (the bias test); measure H1/H3/H6 + map-to-truth
  FSC + boundary ringing + residual histories; the cheap locres-diagonal
  control (§11.3) runs in the same sweep.
- Abort criterion: if omitted-domain damage exceeds its bound at every
  useful strength, stop; the §11.2 aux-competition validator becomes a
  prerequisite, not an option.

### Stage 5 — Gate C science, real data + cost

- **Primary real-scenario harness: abinitio3D with `automsk=yes` in the late
  stages, prior flags on vs off, everything else identical.** Because the
  toggles ride the existing refine3D plumbing (R10), the A/B is two command
  lines differing only in prior strengths — a direct comparison in the exact
  configuration applications will use. Read: final map quality and
  resolution, the DIAG trajectory, the H6 low-shell/centre diagnostics from
  a post-run `rec3D_backends`, and the solvent RMS diagnostics.
- Standalone sweeps on bgal and streptavidin at converged settings against
  matched strength-zero runs (fixed poses, isolates the estimator); then H4
  at the production budget, including the warm-start interplay (the prior
  changes the map the next iteration warm starts from).
- Reliability control (R7): a 10x abinitio3D batch at the chosen strength vs
  the recorded intrinsic ~1/10 rate; treat small-n differences as noise
  unless they replicate.

### Stage 6 — Wilson spectrum source

Sequential after solvent (one variable, R1), reusing the `P_tau` mechanism
(§5.2):

- Gate A addition: swapping spectrum *source* with an identical spectrum is
  bit-identical (mechanism/spectrum separation proven); PSD and curvature of
  the summed prior action; synthetic nonzero-mean RHS correctness.
- Gate B as stage 3; then the four-way ablation (A/B/C/D) at converged
  settings with fixed poses, envelope, spectrum, and band; Gate D criteria
  from §8 decide experimental-release status.

### Stage 7 — explicitly out of first scope

Each item gated behind a successful D run: the NU aux-competitor validator
slot (§11.2, promoted earlier only by the Stage-4 abort); the locres
nonstationary spectral prior (§11.3, judged against post-hoc NU filtering);
`dmats` graded confidence (§11.4); preconditioning of `lambda_s D_s`
(rank-one term stays exact; own parity study); Singer/Gilles off-diagonal
covariance; the lagged LocScale surrogate.

## 11. The NU machinery as the prior infrastructure

The nonuniform-regularization machinery
(`doc/policies/nonuniform_filtering_policy.md`, `src/main/nu_filt/`) is the
mature evidence engine this proposal should build on — the LocScale analogy
made concrete, with better provenance: LocScale-2.0 derives local confidence
from a pseudoatomic reference (model bias), while NU derives it from
cross-validated half-map evidence with a null model. All feeds below are
lag-one by construction (NU runs after the solves inside `volassemble`), the
same discipline as the envelope. Ranked:

1. **Envelope -> solvent weights** (this proposal, §4): the established "in".
2. **Aux-competition -> per-voxel prior validator.** Feed the *priored* ML
   pair through the auxiliary-replacement mechanism (as a distinct
   "competitor" slot, respecting the resolution gate and the `nu_refine`
   exclusion) so the cross-validated NU selection decides per voxel whether
   the priored reconstruction beats the base. The prior then survives only
   where half-map evidence supports it — converting the "envelope misses weak
   density" risk from a hand-designed mask-perturbation test into an
   automatic runtime guard built from trusted machinery. Recommended as a
   design principle for Gate C/D.
3. **Local-resolution field -> nonstationary spectral prior.** The
   Potts-smoothed label field `k_c(v)` (`_nu_locres`) defines the §5
   local-covariance prior with NU label regions in place of LocScale
   windows: `Q_locres = sum_j M_j^T Hp_j^T Hp_j M_j` (soft label-region
   selectors, high-pass above each region's cutoff; every summand PSD, fixed
   during the solve). Cost ~one FFT pair per retained label per matvec
   (8-16 labels); the honest question is the increment over post-hoc NU
   filtering — in-solve, the prior participates in the deconvolution and
   conditions the ML replay rather than masking output — and the A-vs-D
   ablation answers it. A cheap degenerate control belongs in the sweep:
   locres-modulated real-space diagonal Tikhonov (damp voxels with coarse
   evidenced resolution), free per matvec.
4. **Unary Huber evidence (`dmats`) -> graded confidence.** Second-generation
   solvent weights `w^2 = p^2 (1 - c(v))` generalize flatness to "proportional
   to lack of evidence"; the pull-to-common-mean coupling stays restricted to
   genuine solvent (weak in-particle voxels get variance damping without the
   mean term). After the binary envelope form validates.
5. **Potts machinery -> smooth every prior field.** Any per-voxel weight the
   solver consumes should pass through the existing beta-estimated
   ordered-label smoothing so prior fields are spatially coherent. Direct
   reuse, nearly free.
6. **Later:** the `nu_refine` frontier-challenge ratchet as evidence for how
   hard `P_tau` clamps individual beyond-band shells; and the NU §6
   null-estimation discipline reused verbatim for the prior's solvent
   mean/RMS diagnostics.

Caution: with envelope, locres, and evidence all derived from priored maps
and fed back, the §9 support-shrinkage risk generalizes to
**resolution-field drift**. Extend occupancy/overlap tracking to the locres
field, and extend the omitted-domain test to confirm a domain absent from the
locres field can still recover.

## 12. Recommended reading order

Wilson 1942 [3]; Wilson 1949 [4]; Singer 2021 [5]; Gilles and Singer 2022
[6]; Wang 1985 [7]; Terwilliger 1999 [9]; Terwilliger et al. 2020 [10];
LocScale-2.0 [1].

## References

1. A. Bharadwaj, R. de Bruin, and A. J. Jakobi, "Confidence-guided cryo-EM
   map optimisation with LocScale-2.0," *Nature Communications*, vol. 17,
   article 8778, 2026. [doi:10.1038/s41467-026-75327-8](https://doi.org/10.1038/s41467-026-75327-8)

2. A. J. Jakobi, M. Wilmanns, and C. Sachse, "Model-based local density
   sharpening of cryo-EM maps," *eLife*, vol. 6, e27131, 2017.
   [doi:10.7554/eLife.27131](https://doi.org/10.7554/eLife.27131)

3. A. J. C. Wilson, "Determination of absolute from relative X-ray intensity
   data," *Nature*, vol. 150, p. 152, 1942.
   [doi:10.1038/150152a0](https://doi.org/10.1038/150152a0)

4. A. J. C. Wilson, "The probability distribution of X-ray intensities,"
   *Acta Crystallographica*, vol. 2, pp. 318-321, 1949.
   [doi:10.1107/S0365110X49000813](https://doi.org/10.1107/S0365110X49000813)

5. A. Singer, "Wilson statistics: derivation, generalization and applications
   to electron cryomicroscopy," *Acta Crystallographica Section A*, vol. 77,
   pp. 472-479, 2021.
   [doi:10.1107/S205327332100752X](https://doi.org/10.1107/S205327332100752X)

6. M. A. Gilles and A. Singer, "A molecular prior distribution for Bayesian
   inference based on Wilson statistics," *Computer Methods and Programs in
   Biomedicine*, vol. 221, article 106830, 2022.
   [doi:10.1016/j.cmpb.2022.106830](https://doi.org/10.1016/j.cmpb.2022.106830)

7. B.-C. Wang, "Resolution of phase ambiguity in macromolecular
   crystallography," *Methods in Enzymology*, vol. 115, pp. 90-112, 1985.
   [doi:10.1016/0076-6879(85)15009-3](https://doi.org/10.1016/0076-6879(85)15009-3)

8. K. D. Cowtan and P. Main, "Improvement of macromolecular electron-density
   maps by the simultaneous application of real and reciprocal space
   constraints," *Acta Crystallographica Section D*, vol. 49, pp. 148-157,
   1993. [doi:10.1107/S0907444992007698](https://doi.org/10.1107/S0907444992007698)

9. T. C. Terwilliger, "Reciprocal-space solvent flattening," *Acta
   Crystallographica Section D*, vol. 55, pp. 1863-1871, 1999.
   [doi:10.1107/S0907444999010033](https://doi.org/10.1107/S0907444999010033)

10. T. C. Terwilliger, O. V. Sobolev, P. V. Afonine, P. D. Adams, and
    R. J. Read, "Density modification of cryo-EM maps," *Acta
    Crystallographica Section D*, vol. 76, pp. 912-925, 2020.
    [doi:10.1107/S205979832001061X](https://doi.org/10.1107/S205979832001061X)

11. D. Kimanius, G. Zickert, T. Nakane, J. Adler, S. Lunz, C.-B. Schönlieb,
    O. Öktem, and S. H. W. Scheres, "Exploiting prior knowledge about
    biological macromolecules in cryo-EM structure determination," *IUCrJ*,
    vol. 8, pp. 60-75, 2021.
    [doi:10.1107/S2052252520014384](https://doi.org/10.1107/S2052252520014384)
