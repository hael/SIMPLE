# PCG priors: experiment record and the NU-evidence prior design

## Status

This is the single reference for regularization of the PCG reconstruction
backend: what exists on master, what was tried and removed (with the
measurements that killed it), and the approved development order for a direct
NU-evidence prior followed by a Wilson molecular prior. It consolidates and
replaces `pcg_euclid_crash_investigation.md`,
`pcg_reconstruction_production_readiness.md`,
`pcg_refine3D_integration_plan.md`, and
`pcg_real_space_solvent_flatness_prior_proposal.md`. Current implemented
contracts live in `doc/policies/reconstruct3D_pcg_policy.md` and the PCG
section of `doc/policies/refine3D_policy.md`.

**Decision (2026-08-27): the direct NU-evidence prior is priority 1.** When NU
regularization is active it replaces the FSC/SSNR `P_tau` replay precision; it
is never added to `P_tau`. The binary NU-envelope solvent prior is scientifically
retired and must be removed from the target workflow. Its implementation and
measurements remain recorded below only as experiment history until the code is
removed. Wilson is priority 2, after the NU prior has passed its gates.

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
- **Implemented solvent experiment, pending removal**: in NU-enabled PCG runs
  with a valid `nu_envmask3D_stateNN.mrc`, master can also add the binary-
  envelope solvent precision `lambda_s Q_s` to the ML replay. The calibrated
  control currently defaults to `pcg_solvent_lambda_rel=0.1`. This is not the
  approved target architecture after the 2026-08-27 decision; §4 and the
  completed Stage 2--5 entries retain its algebra and measurements as history.
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

Prior priority (decision 2026-08-27):

| Idea | Priority | Recommendation |
| --- | --- | --- |
| direct NU-evidence-conditioned precision | **1 — immediate** | replace `P_tau` in NU mode; §3 and §5 |
| Wilson molecular precision | **2 — after NU** | separate follow-on experiment; §5.5 |
| binary-envelope solvent precision | retired | remove from the target workflow; preserve §4 as experiment history |
| soft state weights | integration semantics | weight both `B` and `D` |
| symmetry projection/permutation | performance only | after distributed profiling |
| total variation / non-negativity / wavelet L1 | deferred | separate nonlinear/proximal solver project |
| joint pose/shift refinement | deferred | separate program (`continuous_3D_refinement_on_pcg_operator.md`) |

Any prior must preserve half independence, positive-semidefinite CG structure,
conditioning, data-mass scaling, and separate reporting of fit and prior
energy, and must be evaluated against a strength-zero control.

## 3. Approved architecture: mode-exclusive regularized replay

The regularized replay has one precision source selected by filtering mode.
There is no additive stack of global ML, solvent, Wilson, and NU priors.

### 3.1 Ordinary PCG mode

Outside NU mode, keep the production estimator unchanged:

```text
base:    (H_data + lambda_0 I) x = b
replay:  (H_data + P_tau + lambda_0 I) x = b
```

`P_tau` remains the global FSC/SSNR shell precision derived from the current
unregularized half pair and raw `D`. Nothing in the NU experiment changes or
reinterprets this path.

### 3.2 NU mode

When the NU machinery is active, replace `P_tau` with one NU-evidence-
conditioned precision `Q_NU(E)`:

```text
base:    (H_data + lambda_0 I) x = b
evidence E = NU analysis(base_even, base_odd)
replay:  (H_data + Q_NU(E) + lambda_0 I) x = b
```

The sequence is current-iteration empirical Bayes, exactly as the existing ML
replay derives `P_tau` from the current base half pair. The NU evidence state is
built after both base solves, then frozen before either replay. It contains no
phase-bearing target and changes only the precision, never the RHS.

Fixed rules:

1. Derive `E` only from the unregularized `_unfil` even/odd pair. Do not let an
   ML-priored auxiliary replacement candidate validate or parameterize its own
   prior.
2. Keep the base maps, FSC/cFAR, raw `(B,D)`, and resolution authority
   unchanged.
3. Attach exactly one replay precision. `P_tau` and `Q_NU` are mutually
   exclusive; a hard assertion must reject simultaneous attachment.
4. Generate no binary molecular envelope for PCG and take no complement to
   manufacture solvent confidence. There is no `Q_s` in the NU target path.
5. Preserve the spherical `mskdiam` support needed by the NU objective and its
   noise/null estimation. Removing envelope generation does not make the NU
   evidence domain self-defining.
6. Keep the NU selected-cutoff/local-resolution products as evidence
   diagnostics and, where required, LP-set bandwidth handoffs. Initially write
   post-hoc `_nu_filt` maps only as controls; the in-solve NU replay is the
   candidate matching reference.
7. Preserve one base solve plus one replay per half and reuse the same raw
   accumulation. No particle reread or third solve is introduced.
8. Wilson is not mixed into this first experiment. It begins only after the
   direct NU estimator passes Gate D.

Two verification prerequisites apply to the NU operator:

- **The `P H P` contract.** The algebra below assumes the support operator is
  applied on both sides of the normal operator and once on the RHS
  (`reconstruct3D_pcg_policy.md` §4 documents this as the contract); verify
  the shared `apply_normal` wrapper and RHS construction actually satisfy it
  before building on it.
- **The deapodization domain.** Any spatially varying NU precision is defined
  on the deapodized map domain. Restrict it to deapodized kernel mode (hard
  error otherwise), as the solvent experiment and `pcg_lambda_lap` did.

## 4. Retired binary-envelope solvent experiment

This section is an experiment record, not an implementation proposal. The
operator was implemented and calibrated successfully, but the 2026-08-27
decision supersedes its binary molecular/solvent partition with the direct NU
precision of §5. Retain the algebra, gates, and real-data measurements so the
experiment is not repeated; remove `Q_s`, its envelope loader, controls, and
convergence guidance from the target workflow when the NU replacement lands.

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
  an expected or reference spectrum (not PCG-compatible as-is; Stage 8).
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

A historical solvent-experiment hypothesis was that the known PCG low-shell
anomalies — the k=1–2 backend amplitude excess (~1.3–2.0) and the centre-bin
radial deficit — are plausibly solvent-offset/low-k noise behavior that `Q_s`
addresses directly. The gated `rec3D_backends` shell and centre-bin
diagnostics measure this for free. Re-baseline those numbers after the ML
warm start is validated, since it moves them too.

## 5. Direct NU-evidence-conditioned precision — priority 1

### 5.1 Why the full NU state is the prior input

The NU engine already evaluates a noise-whitened cross-half prediction cost
`C_c(v)` for every retained local low-pass candidate. That full curve contains
two kinds of information required by a reconstruction prior:

- whether any reproducible signal is supported at voxel `v`;
- the finest scale to which that signal remains supported, including the
  confidence/ambiguity of the scale assignment.

The current envelope path throws most of this information away. It reduces the
bank to `C_20A - min_c C_c`, thresholds the scalar against a whole-support
median/MAD null, solves a binary MRF, filters components, grows the result, adds
a soft edge, and finally takes a complement for `Q_s`. The direct prior must
consume a compact representation of the full candidate curve before that
binary collapse.

The existing bank still lacks one state needed for this use: the coarsest
20-A label absorbs both solvent and genuinely coarse ordered density. Add an
explicit **no-reproducible-signal** candidate or equivalent calibrated null
score. The resulting latent state is

```text
z(v) in {null, 20 A, 15 A, 12 A, 10 A, 8 A, 6 A, 5 A, 4 A, accepted extensions}
```

Its null-versus-signal evidence must remain cross-half predictive rather than a
plain density threshold. Candidate likelihoods/confidences must be calibrated
against the exact bank, whitening profile, and smoothing scales; raw Huber
costs are not automatically normalized posterior probabilities.

### 5.2 Compact evidence state

Do not persist `dmats_mask(n_vox,n_candidates)` between reconstruction phases.
Before releasing it, reduce it to the sufficient state needed by the solver:

- `a_b(v)`: monotone support/confidence that detail band `b` is reproducible;
- the selected local cutoff for diagnostics and LP-set handoff;
- an uncertainty measure such as best-versus-next margin or label entropy;
- the null calibration/provenance needed to reproduce the mapping.

Use the existing ordered-label spatial model to regularize this state, but do
not force it through binary topology, component removal, dilation, or a mask
edge. Confidence must remain graded through molecular boundaries and ambiguous
domains.

### 5.3 Quadratic precision

Let `B_b` be a normalized detail/band-pass analysis operator and let `W_b` be a
nonnegative spatial weight derived monotonically from lack of evidence for
band `b`. Define

```text
R_NU(x | E) = 1/2 sum_b || sqrt(W_b) B_b x ||^2
Q_NU(E)     = sum_b B_b^T W_b B_b
```

Every summand is positive semidefinite, so the fixed replay operator remains
symmetric positive definite after `lambda_0 I` and the existing `P H P`
support sandwich. The detail bank should exclude the global DC mode, making a
constant field a null mode of the NU precision without constructing a solvent
mean explicitly.

The one field has the required limiting behavior:

- null evidence: penalize all non-DC detail, which yields solvent flattening;
- coarse ordered density: preserve supported coarse bands and suppress finer
  unsupported bands;
- strong ordered density: preserve detail through the evidenced local cutoff;
- uncertain boundaries: grade the precision continuously instead of choosing
  a molecular/solvent side.

The initial implementation should use three or four broad, normalized detail
bands rather than one operator per NU label. Record the frame normalization so
adding bands does not silently increase total strength. A full 8--16-label
operator is justified only after the compact form demonstrates scientific gain
and acceptable replay cost. The preconditioner needs a declared nonnegative
approximation to the NU precision; scientific comparisons must still be made
at matched convergence, not matched iteration count.

### 5.4 Relationship to global ML regularization

`Q_NU` is the nonstationary empirical-Bayes alternative to `P_tau`, not a
correction layered on top. The existing global ML path maps shell FSC to SSNR
and scales its precision by shell-mean raw `D`. The NU development must define
and validate the corresponding mapping from calibrated local evidence to
bandwise precision. Until that mapping is established, the evidence field is
a detection/selection statistic, not yet a quantitatively normalized prior.

The mode assertion is load-bearing:

```text
ordinary mode: P_tau present, Q_NU absent
NU mode:       P_tau absent,  Q_NU present
```

### 5.5 Wilson molecular prior — priority 2

Wilson supplies population information about molecular Fourier covariance,
especially where the experiment is uninformative. NU evidence cannot contain
information the data never observed, but the first question is whether such
extrapolation is needed at all: the conservative NU estimator regularizes only
the degrees of freedom the half-map evidence does not support.

Only after the direct NU estimator passes Gate D should a Wilson experiment
begin. Keep it separate from the NU acceptance experiment. The first Wilson
form remains a zero-mean, shell-diagonal precision
`Q_W = F* diag(q_W(k)) F`, implemented through one declared spectrum-source
path rather than a parallel Fourier array. Compare Wilson against the accepted
NU estimator; do not combine them until each has an independently established
benefit and a separate combination experiment is approved.

Singer derives the non-diagonal covariance of Fourier coefficients from the
random bag-of-atoms model [5]; Gilles and Singer use the corresponding mean and
covariance as a Bayesian molecular prior [6]. Off-diagonal structure, nonzero
means, and LocScale-informed local covariances remain later experiments. A
nonzero or phase-bearing mean is outside the first Wilson scope.

## 6. Workflow, evidence contract, and bias discipline

### 6.1 Per-state sequence (one particle accumulation, two solve phases)

1. Reduce the raw, data-only even and odd statistics as today.
2. Base-solve each half unchanged, with neither `P_tau` nor `Q_NU`; retain and
   write the `_unfil` pair.
3. Compute FSC/cFAR from that pair. It remains the unregularized resolution
   authority and is not used as a second NU replay precision.
4. If NU mode is inactive, build `P_tau` and run the existing ML replay.
5. If NU mode is active, run the NU candidate/null analysis on the current
   base pair, compact and freeze the evidence state, build `Q_NU`, assert that
   `P_tau` is absent, and replay both halves with the same fixed precision.
6. Write the standard replay maps. Retain selected-cutoff/local-resolution and
   post-hoc NU products as diagnostics during validation; generate no PCG
   solvent envelope.

### 6.2 Evidence contract

The direct prior accepts only evidence derived from the current state's
unregularized even/odd base pair. Validate before replay:

- identical dimensions, sampling, and physical extent for the base halves;
- finite, nonnegative candidate costs and a valid radial whitening profile;
- a nonempty spherical `mskdiam` support with an adequate null population;
- exact candidate-bank, smoothing, extension, and null-calibration provenance;
- monotone band-support confidences in `[0,1]` and finite nonnegative precision
  weights;
- one immutable evidence identity used for both half replays.

Emit at least:

```text
pcg_replay_prior_mode=global_ml|nu_evidence
pcg_nu_prior_enabled=  pcg_nu_candidate_count=  pcg_nu_null_fraction=
pcg_nu_supported_fraction_bandNN=  pcg_nu_uncertain_fraction=
pcg_nu_prior_energy_final=  pcg_nu_preconditioner_mode=
pcg_nu_evidence_source=base_unfil  pcg_nu_evidence_provenance=
```

There is no mask loader, missing-envelope skip, mask resampling, morphology,
or silent fallback in the NU prior path. Failure to construct valid evidence
must be explicit; fallback policy, if any, is a workflow decision rather than
an inferred substitute prior.

### 6.3 Half-set and bias discipline

Both halves use the same fixed scalar confidence/bandwidth field, so their
errors can still become correlated even though the prior has no phase-bearing
mean. The reported FSC must therefore remain that of the unregularized base
pair. Keep raw accumulators half-specific and unchanged; freeze the evidence
before either replay; report shipped-pair FSC only as a correlation-inflation
diagnostic; compare against the ordinary `P_tau` estimator using map-to-truth,
independent map-to-model, or held-out predictive evidence; and track local-
cutoff/confidence overlap across refinement iterations. Deriving the field
from the current base pair removes lagged envelope feedback, but matching
against the NU-regularized replay map remains an outer refinement feedback
loop and weak omitted domains must be tested explicitly.

### 6.4 Ownership

| Concern | Owner |
| --- | --- |
| Mode selection/defaults/validation | `parameters`, dictionaries, UI, and refinement policy; ordinary global ML versus NU replay is explicit and mutually exclusive |
| Numerical precision | `src/main/volume/simple_reconstructor_pcg.f90` (fixed NU evidence state behind the shared normal-operator wrapper and inside the support sandwich) |
| Evidence construction | `src/main/nu_filt/` (null candidate, full unary reduction, spatial regularization, compact band-confidence state) |
| PCG reconstruction orchestration | `simple_rec3D_pcg_strategy.f90` and assembly owner (base pair -> FSC -> evidence -> one replay; shared/distributed parity) |
| NU diagnostics/LP handoff | `simple_commanders_rec_distr.f90` and NU modules; avoid recomputing an inconsistent second evidence state |
| Policy record | `reconstruct3D_pcg_policy.md`, `nonuniform_filtering_policy.md` |

The raw `(B,D)` format does not change. NU evidence and precision are
solve-time state, never accumulated into `D`; constructing the prior triggers
no particle or raw-artifact reread.

## 7. Hypotheses

- **H1:** A single NU-conditioned replay suppresses unsupported solvent and
  high-resolution variation without a binary molecular/solvent partition.
- **H2 (architectural invariant):** Selecting NU replay changes neither raw
  `(B,D)`, the base solves, nor FSC/cFAR, and preserves one base solve plus one
  replay per half.
- **H3:** The explicit null state preserves genuinely coarse ordered density
  better than treating the 20-A saturating label as solvent.
- **H4:** `Q_NU` produces at least the scientific benefit of post-hoc NU
  filtering while participating in deconvolution and improving or preserving
  replay conditioning.
- **H5:** NU replay outperforms the ordinary global `P_tau` replay on
  heterogeneous/local-resolution data without degrading uniform high-SNR
  cases materially.
- **H6 (performance gate):** The compact band operator and its preconditioner
  reach matched residual accuracy within an operationally acceptable budget.

Record numerical thresholds for "material loss", parity, convergence, and
acceptable overhead *before* the real-data sweep.

## 8. Staged validation

### Gate A: algebra and operator identity

Small synthetic volumes and fixed evidence fields: linearity of `Q_NU`; the
dot-product identity `<x,Q_NU y>=<Q_NU x,y>`; nonnegative quadratic form;
constant/DC null behavior; zero action for a fully supported band; monotone
action as evidence is withdrawn; continuity under small confidence changes;
penalty gradient versus finite differences; correct `P(H+Q_NU)P`
composition; hard failure when `P_tau` and `Q_NU` are both requested; and
unchanged ordinary-mode numerical identity. Mutation tests must break the
adjoint factorization, band ordering, mode exclusion, and evidence freeze.

### Gate B: workflow and artifact invariants

Selecting NU replay does not change raw `(B,D)` checksums; leaves the base
solves and FSC bit-identical; performs exactly one replay per populated half;
derives evidence from `_unfil` maps without an ML auxiliary candidate; freezes
one evidence identity for even and odd; never attaches `P_tau`; generates no
PCG solvent envelope; keeps shared/distributed parity under the existing
deterministic-reduction tolerance; and leaves ordinary global-ML mode
bit-identical.
(Fractional/trailing reconstruction is hard-errored on the PCG path today,
which keeps this gate's surface small; revisit when that guard lifts.)

### Gate C: controlled scientific tests — with convergence isolation

The direct lesson of the removed nu^4 prior is that, at production iteration
budgets, measured "prior effects" are confounded with replay under-convergence.
The spatially varying NU precision is not representable exactly by the current
Fourier-diagonal preconditioner, so convergence cost must be measured.
**Scientific comparisons run at converged settings** (high `maxits_pcg` /
tight `rtol`); the production-budget behavior is measured separately as H4.

Harness: the neutral-phantom fixture and the gated `rec3D_backends`
ground-truth mode (map-to-truth FSC, radial LS profiles, background handling,
centre-bin and shell diagnostics) already provide the measurement
infrastructure. Measure map-to-truth FSC (not only half-map FSC), local
resolution calibration, background variation, boundary ringing/leakage,
recovery of deliberately weak and coarse peripheral domains, PCG residual
history/stop reason, alignment overlap, and shipped-versus-base half-map
correlation inflation.

The first ablation is deliberately two-way:

| Run | Replay precision | Purpose |
| --- | --- | --- |
| A | existing global FSC/SSNR `P_tau` | Production estimator/control. |
| B | direct `Q_NU`, with `P_tau` absent | Isolate the complete NU estimator. |

Use identical particles, poses, raw accumulators, base maps, candidate bank,
and converged solver settings. Within B, sweep only one declared NU precision
temperature/strength if calibration has not eliminated it. Do not tune from
half-map FSC alone; prefer held-out predictive likelihood, map-to-truth FSC,
or independent map-to-model comparison. Wilson is not part of this gate.

### Gate D: acceptance for an experimental release

All algebraic and workflow invariants pass; ordinary global-ML outputs are
unchanged; improvements are supported by independent truth/model or held-out
evidence; coarse and weak-domain recovery is not systematically worse;
background suppression does not require a binary envelope; the NU replay
outperforms or sits on a useful quality/cost frontier versus `P_tau`; per-half
convergence stays within the operational budget or a measured new budget is
approved; and shared/distributed routes agree. Only then does Wilson become
the active development priority.

## 9. Risks and mitigations

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| Null state fails to separate solvent from coarse density | Real coarse signal is over-regularized. | Explicit no-signal competitor; coarse-domain synthetic gate; calibrate against the exact bank. |
| Huber costs treated as posterior probabilities | Arbitrary and nonportable precision scale. | Calibrated null/temperature with recorded bank, whitening, and smoothing provenance. |
| Evidence-derived replay map drives later alignment | Self-reinforcing local-bandwidth loss. | Current base-pair evidence, full spherical support, weak/omitted-domain recovery tests, field-overlap tracking. |
| Common evidence precision inflates half-map correlation | Misleading shipped-pair FSC. | Resolution from base pair only; shipped FSC diagnostic only; independent truth/model evidence. |
| `P_tau` accidentally remains active | Uninterpretable combined estimator. | Hard mutual-exclusion assertion and mutation test. |
| Wilson is introduced before NU is resolved | Multiple moving scientific variables. | Wilson starts only after NU Gate D. |
| LocScale amplitude target treated as linear | Operator depends on `x`; PCG assumptions fail. | Fixed Wilson covariance or a lagged surrogate; never rescale inside `apply_normal`. |
| Common LocScale target carries phases | Half errors directly correlated. | Share only phase-free profiles; lagged, half-specific otherwise. |
| Overlapping detail bands change total strength | Candidate-bank size silently changes regularization. | Declared frame/band normalization; parity test when refining the bank. |
| Prior effect confounded with under-convergence | Wrong scientific conclusions (the nu^4 lesson). | Gate C convergence isolation; H4 measured separately. |
| Preconditioner cannot represent spatial `Q_NU` exactly | Iteration-count regression at production budgets. | Nonnegative declared approximation; matched-convergence science; cost gate. |
| A second NU pass builds different evidence | Replay and postprocess disagree about local support. | Construct once from base halves; share the compact state with diagnostics/LP handoff. |

## 10. Staged development plan

Each active stage ends at a gate; a later stage does not begin until the
preceding gate passes and the stage's changes are committed. New NU behavior
must be explicit and protected while the existing ordinary global-ML path
remains unchanged. "Fixture" means the neutral-phantom project with ground
truth plus the gated `test=rec3D_backends`. Stages 2--5 below are retained as
the completed/superseded solvent experiment record, not as active work.

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
- **R8 — evidence plumbing before operator math.** First expose the null
  candidate and compact immutable evidence state, then attach `Q_NU`.
- **R9 — record acceptance thresholds in this document before running the
  experiment that is judged by them.**
- **R10 — replay precisions are mode-exclusive.** Ordinary PCG attaches
  `P_tau`; NU PCG attaches `Q_NU`. No command-line combination may activate
  both. Wilson is a later estimator choice, not an additive knob in the NU
  acceptance experiment.

### Stage 0 — historical warm-start foundation

The cross-iteration ML warm start changed the baselines against which the
solvent experiment was measured (R2). These entries are retained as history;
they do not block the active Stage 6 NU design.

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
  `report_solvent_envelope_status` in the strategy runs the former mask contract
  (presence, cubic lattice, physical-extent identity, constant-FOV
  `read_and_crop` resample on lattice mismatch, finiteness, `[0,1]` re-clip,
  nonzero solvent evidence) and emits the `pcg_solvent_*` block with
  `prior_enabled=F`, once per state before either half's ML replay, both
  execution paths. Fixture gate PASS unchanged (median ratio .8428, identical
  to the pre-change runs). NOTE: the fixture runs `objfun=cc`, and
  `l_ml_reg` requires `objfun=euclid`, so the diagnostic block cannot fire
  there — the envelope-absent skip and the resample demo must be read from
  the next production refine3D/abinitio3D run.

### Stage 2 — historical solvent `Q_s` operator record

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

### Stage 3 — historical solvent workflow-integration record

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
- **Harness-validated (2026-08-25, neutral fixture, see the alignment-free
  loop below):** prior-on at strength 1e-3 — `prior_enabled=T`,
  `skip_reason=none`, per-half `lambda_eff`/`mean`/`rms`/`penalty` lines,
  gates PASS (in-band .9904 vs .9927 unpriored; centre-bin ratio 1.16);
  strength-zero control — `skip_reason=strength_zero`, .9927 identical to
  the no-envelope run (R5); envelope-absent at positive strength —
  THROW_WARN with `skip_reason=envelope_absent`, run completes and PASSES.
- **Pending on real data** (user-side run): one full abinitio3D with
  `rec_backend=pcg automsk=yes` and a positive strength completing, showing
  the prior activating exactly when the envelope first exists lag-one and
  exercising the stage-handoff resample.

### DECISION (2026-08-26): automasking auto-enabled for pcg in abinitio3D

With the prior on by default, abinitio3D closes the loop itself: when
`rec_backend=pcg` and NU filtering is active (the default), the stage policy
forces `automsk=yes` from the first NU stage (`NU_FILTER_STAGE`) onward — the
first NU stage's assembly generates the state-specific evidence envelope and
every later stage's ML replays are solvent-priored lag-one, with no user
flags. An explicit `automsk=no` vetoes; a user-supplied mode (`tight`) is
respected; gridding runs keep the last-stage-only automasking behavior.
Verification (pending, user-side): one abinitio3D run with `rec_backend=pcg`
showing NU ENVELOPE OCCUPANCY at stage 6 and `pcg_solvent_prior_enabled=T`
from stage 7 — this now doubles as the outstanding Gate B lag-one activation
item.

### DECISION (2026-08-26): solvent suppression readout in the convergence report

With the prior default-on and self-bootstrapping, refinement needs an in-run
answer to "is the prior firing, and is the strength right?" that requires no
harness and no ground truth. The readout is the **solvent suppression
percentage**, per half:

    supp% = 100 * (1 - rms_s(x_ML) / rms_s(shrink(x_base)))

where `shrink(x_base)` is the base solution passed through the closed-form
per-shell FSC/Wiener shrinkage (`fsc2shrink_filter` in `simple_estimate_ssnr`,
shared with the regularized ML initialization). The shrunk base is the
reference rather than the raw base because the Wiener component of the ML
replay reduces solvent variation on its own; referencing it out isolates the
solvent prior's contribution, so an inert prior reads ~0% regardless of tau.
Plumbing: `report_solvent_solve_stats` prints `pcg_solvent_rms_ref` and
`pcg_solvent_suppression_pct` per half; both PCG execution paths average
even/odd per state and persist `simple_pcg_solvent_stats.txt`
(delete-then-rewrite each ML volassemble, so a skipped prior leaves no stale
values); `check_conv3D` prints `% SOLVENT SUPPRESSED (PCG SOLVENT PRIOR)`
with the other iteration stats, writes it to `simple_stats.txt`
(`PCG_SOLVENT_SUPPRESSION_PCT` + per-state keys), and advises:

- **< 5%**: prior inert — increase `pcg_solvent_lambda_rel` (~3x)
- **5-60%**: nominal — keep the strength
- **> 60%**: over-flattening risk — decrease (~3x)

Thresholds are provisional constants in `simple_convergence`
(`PCG_SOLVENT_SUPP_INERT/OVER_PCT`), anchored to the bgal ladder verdicts;
recalibrate by rereading the ladder now that every run prints the readout.
**Automatic sweep (2026-08-26):** `test=rec3D_backends` now sweeps the
default ladder {0, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1} automatically whenever a
prior-capable invocation (ml_reg + NU envelope in the cwd) leaves
`pcg_solvent_lambda_rel` unset -- one execution directory per rung, collated
Stage-5 observables printed and written to `rec3D_backends_sweep_summary.txt`
(base/shipped pair FSC, suppression %, solvent RMS, gated-band ratio/FSC,
radial min/max, truth FSC, gates), with the base-pair-spread negative control
asserted across the ladder. An explicit strength runs a single comparison --
the follow-up mechanism for extra rungs.

Fixture sweep record (its=5, rtol=1e-3, lp=10, reg-init; base pair pinned
4.036/3.347 across the ladder, spread 0.0000):

| lambda | supp% | sol_rms | ship 0.5/0.143 | band_ratio | band_fsc | rad_min/max | truth_fsc_pcg |
|---|---|---|---|---|---|---|---|
| 0     | --    | --   | 4.036/3.347 | 1.006 | 1.000 | 0.98/1.01 | 0.971 |
| 1e-3  | 5.0   | 12.2 | 4.036/3.347 | 1.002 | 1.000 | 0.99/1.01 | 0.970 |
| 1e-2  | 14.7  | 11.  | 3.920/3.347 | 0.987 | 0.999 | 0.97/1.06 | 0.970 |
| 3e-2  | 29.1  | 8.8  | 3.920/3.347 | 0.950 | 0.997 | 0.92/1.12 | 0.966 |
| 1e-1  | 51.7  | 6.0  | 3.709/3.267 | 0.872 | 0.981 | 0.86/1.23 | 0.948 |
| 3e-1  | 71.2  | 3.6  | 3.430/2.744 | 0.798 | 0.930 | 0.81/1.42 | 0.898 |
| 1     | 86.1  | 1.7  | 2.144/2.144 | 0.759 | 0.848 | 0.75/1.85 | 0.810 (1 gate) |

Fixture read (R6: judged by truth, not half-map FSC): truth agreement is
flat to 1e-2, mildly down at 1e-1 (-0.023) and clearly damaged from 3e-1;
first material truth loss coincides with supp ~52% and clear damage with
~71%, bracketing the provisional 60% over-flattening bound from both sides.
The fixture's truth-optimal strength (~1e-2) is below bgal's visually-chosen
1e-1 -- dataset-dependent, which is exactly what the in-run suppression
readout is for. The 5% inert bound is confirmed: the 1e-3 rung reads 5.0%
with every other observable at the zero-control value.

### DECISION (2026-08-26): solvent prior on by default at the calibrated strength

Following the bgal calibration the solvent prior graduated from experimental
toggle to default behavior: `pcg_solvent_lambda_rel` now defaults to **0.1**
(the calibrated operating strength), so whenever `rec_backend=pcg` runs an ML
replay and a valid state-specific NU evidence envelope exists, the prior is
applied — no flag required. `0` remains the explicit off-switch, preserving
the A/B capability R10 requires for the remaining Gate C envelope-variant
tests and for calibration on new dataset classes. With default-on, an ABSENT
envelope is the normal lag-one state (early iterations, gridding-era
projects) and produces only the diagnostic block; a PRESENT-but-unusable
envelope (extent mismatch, invalid values) still warns loudly. The former
never-substitute-a-mask rule remained in force for that experiment.

### Stage 5 record: streptavidin sweep (2026-08-27, minimal-invocation harness)

First bare-invocation sweep on real data (10335 streptavidin, d2, box 128,
smpd 1.072, its=5): base pair pinned 3.430/3.119 across the ladder, spread
0.0000, all gates PASS.

| lambda | supp% | shipped 0.5/0.143 | band_fsc | rad_min | verdict |
|---|---|---|---|---|---|
| 0     | --   | 3.347/3.049 | 0.999 | 1.000 | control |
| 1e-3  | 6.6  | 3.347/3.049 | 0.999 | 1.000 | inert |
| 1e-2  | 21.8 | 3.267/3.049 | 0.996 | 1.000 | clean |
| 3e-2  | 36.1 | 3.191/3.049 | 0.991 | 1.000 | clean -- operating point |
| 1e-1  | 51.2 | 3.119/2.800 | 0.981 | 0.998 | shipped inflation fires |
| 3e-1  | 63.5 | 3.049/2.144 | 0.965 | 0.912 | inflation + rim erosion |
| 1     | 78.7 | 2.144/2.144 | 0.932 | 0.741 | over-flattened |

Cross-dataset findings: (a) the suppression-vs-lambda curve nearly overlays
the fixture's (6.6/22/36/51/64/79 vs 5/15/29/52/71/86 per decade) -- the
data-scale anchoring makes the readout portable; (b) the FIRST
over-regularization signature fires at ~51% suppression on both datasets
(fixture truth-FSC loss at 51.7%, streptavidin shipped inflation at 51.2%),
independent observables agreeing that trouble onset is a property of the
suppression coordinate. OPEN DECISION: lower PCG_SOLVENT_SUPP_OVER_PCT from
60 to ~45-50 (at 60 the default 0.1 reads "nominal" on streptavidin while
its inflation diagnostic has fired), and whether the default strength should
move 0.1 -> 3e-2 (clean everywhere measured) with the convergence guidance
steering upward toward the ~50% line.

### Gate B lag-one activation: CLOSED (2026-08-27, streptavidin abinitio3D)

The warm-started `rec_backend=pcg` abinitio3D run verified the full
self-bootstrapping staging on real data with zero user flags: stage 5
automasking off + `envelope_absent` quiet skip; stage 6 first iteration
forces automasking on, NU engages with the spherical reference ("envelope
not available yet" — correct first-iteration state) and the volassemble
generates the envelope (occupancy 7.4%); the next iteration consumes it
lag-one (`prior_enabled=T`, `skip_reason=none`) and every later iteration
re-generates/re-consumes. Steady-state warm-started suppression: ~28-35%,
stable (third regime datum — cold its=5 box-128: 51%; converged bgal
box-256: 14%; warm box-88: ~32% — triply confirming the removal of
absolute-% upper thresholds). The prior-free base-pair 0.143 improved
3.812 -> 3.191 A through stage 6 (crop Nyquist) — the loop demonstrably
improves registration. WATCH ITEMS (not acted on): at prior/NU engagement
the promoted matching lp (6.2 -> 3.5 A), the switch to envelope-masked
references, and the prior all hit the matcher in one iteration ->
orientation overlap collapses 0.70 -> 0.04 and recovers slowly; stage 6
exits unconverged on its budget. Envelope occupancy oscillates 7 -> 23 ->
14-18% (jitter, not runaway shrinkage). Candidate mitigations: stagger the
reference-mask switch one iteration after lp promotion; cap the promoted
lp a margin below the crop Nyquist.

### DECISION (2026-08-27): inflation-based over-flattening guidance (implemented)

The convergence guidance is restructured per the anchor verdict below:

- **Inert**: suppression < 2% (`PCG_SOLVENT_SUPP_INERT_PCT`) -> increase
  `pcg_solvent_lambda_rel` ~3x. Dataset-robust at the low end.
- **Over-flattening**: shipped-pair FSC=0.143 pulling > 5% finer than the
  base-pair 0.143 (`PCG_SOLVENT_INFL_PCT`) -> decrease ~3x. The portable
  signal (fired at 1e-1 on streptavidin, quiet through 3e-1 on bgal,
  matching both ladder verdicts). Absolute suppression thresholds are NOT
  used for the upper bound.
- Suppression % is retained as the monotone within-run trend readout.

Plumbing: the PCG strategy now measures the shipped (regularized) pair's
FSC crossings in-run after the ML replays (`shipped_pair_res`, same soft-
masked procedure as the harness diagnostic; computed only when the prior
fired) and persists `PCG_BASE_FSC0143_STATEXX`/`PCG_SHIP_FSC0143_STATEXX`
in `simple_pcg_solvent_stats.txt`; `check_conv3D` prints
`% SHIPPED-PAIR FSC INFLATION (PCG PRIOR)` (per-state lines for
`nstates>1`), applies the two-test guidance, and records
`PCG_SHIP_INFLATION_PCT` in `simple_stats.txt`. The shipped-pair crossing
shares regularization between halves and is never a resolution claim.

### Stage 5 record: bgal its=30 anchor at 1e-1 (2026-08-27) — portability verdict

Base solves converge (RESID 2.7e-2 @ 30 its); ML replays stop on XTOL at
18/19 its with RESID still 0.12 (CG stalls by step size on the box-256 ML
system — preconditioner item). Suppression at 1e-1: 12.8/14.0% (up from
4.7% at its=5; rms 0.675 -> 0.638; reference rose to ~0.73 because the
30-it base carries more high-frequency noise). VERDICT: convergence
explains ~3x of the bgal deficit but NOT the gap to streptavidin's 51% —
lambda_eff at identical lambda_rel=0.1 is 1.98e4 (bgal) vs 1.53e6 (strep),
i.e. the s_data(D) anchoring does not equalize prior-vs-data balance across
box sizes; bgal's dose-response is genuinely right-shifted (~lambda_rel>=1
for 50%-class suppression). Side findings: converged 1e-1 rim is pristine
(bin-9 ratio 1.067 vs 0.851 at its=5 — better peripheral preservation);
Nyquist-edge shells k=126-128 blow up (amp ratio up to 17x) over 30 ML its
— S6 beyond-band growth, outside the gated band, irrelevant at production
budget. GUIDANCE REDESIGN (proposed, pending decision): absolute
suppression thresholds cannot be universal; the observable that landed at
the right per-dataset operating point on BOTH datasets is shipped-vs-base
pair inflation (fired at 1e-1 strep, 3e-1 bgal). Restructure convergence
guidance as: inert = suppression ~<2% (dataset-robust); over-flattening =
shipped-pair 0.143 materially finer than base-pair 0.143 (portable, both
already computed in-run); suppression % retained as the monotone
within-run trend readout.

### Stage 5 record: bgal harness sweep + suppression-portability caveat (2026-08-27)

Bare-invocation sweep on bgal (d2, box 256, its=5, shared): base pinned
4.132/3.709, spread 0.0000, all PASS; the rms ladder and shipped-inflation
verdicts reproduce the 2026-08-26 manual ladder exactly (rms
0.716/0.712/0.703/0.675/0.61/0.46; inflation first at 3e-1). BUT the
suppression scale does NOT overlay streptavidin/fixture: bgal reads
-1/-0.5/+0.7/4.7/14/35 % per decade vs streptavidin's 6.6/22/36/51/64/79 %.
Diagnosis: convergence state — bgal ML RESID ~0.17-0.19 at its=5 (box 256,
STOP=maxits) vs streptavidin ~0.03 (box 128); an unconverged iterate has not
expressed the prior, so cold-solve suppression under-reports on large boxes.
The negative values at inert strengths expose a smaller per-dataset
reference bias (shrunk base is not exactly the lambda=0 ML solution; +6.6%
strep, -1% bgal at inert rungs). CONSEQUENCES: (a) the ~50%-onset
universality claim below is RETRACTED as an absolute-percentage rule —
suppression is monotone and useful within a ladder, comparable across
datasets only at matched convergence; (b) threshold/default retuning is ON
HOLD pending the its=30 bgal anchor at 1e-1 (decisive: convergence vs
data-scale anchoring) and the warm-started abinitio3D suppression
trajectory (the accumulated steady-state number is the one the convergence
guidance actually sees); (c) the harness should report/flag RESID next to
suppression when the ML solve stops on maxits far from rtol.

### Stage 5 record: bgal solvent-prior strength calibration (2026-08-26)

Full dose-response on beta-gal (box 256, smpd 1.275, d2, ~4.6k ptcls,
`maxits_pcg=5` with the regularized ML initialization; every run's base-pair
FSC pinned at 4.132/3.709 A -- the negative control across the 1000x sweep):

| lambda_rel | solvent RMS | rim ratio (64-72 px) | shipped FSC 0.5/0.143 | verdict |
|---|---|---|---|---|
| 1e-3 | 0.7156 | 0.862 | 4.132/3.709 | inert |
| 1e-2 | 0.7117 | 0.861 | 4.132/3.709 | inert |
| 3e-2 | 0.7032 | 0.859 | 4.132/3.709 | whisper |
| 1e-1 | 0.6752 | 0.851 | 4.132/3.709 | active, clean |
| 3e-1 | 0.6075 | 0.827 | 4.080/3.709 | first inflation + mild rim loss |
| 1    | 0.4616 | 0.741 | 4.030/3.667 | over-flattening |

Findings: (1) the active regime on real data starts ~1e-1, roughly 30x above
the neutral fixture's (noise-dependent, as the aliasing argument predicts);
(2) the two over-regularization diagnostics -- envelope-boundary rim erosion
in the radial table and shipped-pair FSC inflation -- fire TOGETHER at 3e-1
and grow monotonically, while the base-pair FSC never moves; (3) backend
agreement (gridding vs pcg) improves monotonically with lambda through the
whole range (band k=2-103 -> 104, high-shell excess shrinking), evidence the
prior removes error rather than adding structure (independent algorithms,
shared data). **Recommended default `pcg_solvent_lambda_rel=1e-1`; 3e-1 is
the aggressive bound.** In refinement, prefer the conservative default: the
shipped map feeds next-iteration alignment and the lag-one envelope, so rim
erosion is the entry point of the self-reinforcing support-shrinkage feedback
(S9) that the standalone harness cannot observe. Earlier note: the first
its=5 sweep (pre-initialization) was transient-dominated and inconclusive --
strength conclusions require the regularized init or converged solves (R3).

### Alignment-free prior validation harness (no new abinitio3D runs)

Prior validation and strength sweeps need only a project file with a previous
abinitio3D registration (poses in `ptcl3D`), the last iteration's half maps,
and that run's `sigma2_it_*.star` files — no new alignments. `nu_filt3D` is
the envelope generator: with `nu_envmsk=yes` it builds the identical NU
evidence envelope object volassemble publishes
(`envmask3D_from_lmask` on the NU evidence margin), just under a different
name. The loop, run in a scratch dir containing the sigma2 star files (or in
the final refine3D stage dir of the abinitio3D run):

```text
simple_exec prg=nu_filt3D vol1=<prev_odd.mrc> vol2=<prev_even.mrc> \
    smpd=<smpd> mskdiam=<D> nu_envmsk=yes mkdir=no nthr=<N>
cp outvol_nu_envmask.mrc nu_envmask3D_state01.mrc
simple_test_exec test=rec3D_backends projfile=<proj.simple> pgrp=<pgrp> \
    mskdiam=<D> objfun=euclid ml_reg=yes pcg_solvent_lambda_rel=<X> nthr=<N> \
    [vol1=<truth.mrc> lp=<lp>]
```

Notes: `vol1` to nu_filt3D is the ODD half, `vol2` the EVEN half (nu_filt3D
convention); the envelope is written as `<outvol-basename>_nu_envmask.mrc`
(default `outvol_nu_envmask.mrc`); `mkdir=no` keeps it in the cwd; a
native-box envelope against a cropped reconstruction is handled by the
loader's constant-FOV resample (and demonstrates that path); `X=0` is the R5
control and must be bit-identical; with `ml_reg=yes` the truth LS-profile
check is reported as a diagnostic rather than gated (the shipped maps are
ML-regularized; the gate is calibrated on unregularized maps). The
envelope-tuning knobs (`nu_msk_sig/beta/dens`, `amsklp`) let the same loop
drive the former Gate C eroded/dilated envelope variants without touching the
reconstruction inputs. Every historical sweep value of `pcg_solvent_lambda_rel`
reuses the same envelope file — one `nu_filt3D` call per envelope variant,
one `rec3D_backends` call per strength. **Validated on the neutral fixture
(2026-08-25):** a 20-second single-iteration refine3D seeds
`sigma2_it_1.star` and the half pair when a project has poses but no sigma2
files; the loop then activates the prior end-to-end
(`pcg_solvent_prior_enabled=T`, per-half lambda_eff/mean/rms/penalty lines).

### Stage 4 — retired solvent Gate C plan (do not execute)

- Record thresholds first (R9): "material molecular loss", acceptable
  weighted-solvent-RMS reduction, omitted-domain damage bound.
- Fixture sweep at converged settings; envelope variants: true, eroded,
  dilated, omitted-domain (the bias test); measure the former solvent
  hypotheses plus map-to-truth
  FSC + boundary ringing + residual histories; the former cheap locres-
  diagonal control was planned for the same sweep.
- Abort criterion: if omitted-domain damage exceeds its bound at every useful
  strength, stop; the former aux-competition validator would have become a
  prerequisite rather than an option.

### Stage 5 — retired solvent Gate C plan and completed measurements

- **Primary real-scenario harness: abinitio3D with `automsk=yes` in the late
  stages, prior flags on vs off, everything else identical.** Because the
  toggles ride the existing refine3D plumbing (R10), the A/B is two command
  lines differing only in prior strengths — a direct comparison in the exact
  configuration applications will use. Read: final map quality and
  resolution, the DIAG trajectory, the historical low-shell/centre diagnostics from
  a post-run `rec3D_backends`, and the solvent RMS diagnostics.
- Standalone sweeps on bgal and streptavidin at converged settings against
  matched strength-zero runs (fixed poses, isolates the estimator); then H4
  at the production budget, including the warm-start interplay (the prior
  changes the map the next iteration warm starts from).
- Reliability control (R7): a 10x abinitio3D batch at the chosen strength vs
  the recorded intrinsic ~1/10 rate; treat small-n differences as noise
  unless they replicate.

### Stage 6 — direct NU-evidence replay — active priority 1

6.1 Evidence contract, no solver behavior change:

- add and validate the no-reproducible-signal/null competitor;
- derive compact monotone band-support confidence, cutoff, uncertainty, and
  provenance from the full unary bank before it is released;
- prove that the state is derived only from the `_unfil` pair and is identical
  for both replay halves;
- expose diagnostics and mutation-test candidate ordering/null calibration.

6.2 Operator Gate A:

- implement the normalized three- or four-band `Q_NU` factorization;
- prove adjoint identity, PSD, DC null, monotonicity, finite-difference
  gradient, and `P(H+Q_NU)P` composition;
- add the hard `P_tau`/`Q_NU` mutual-exclusion assertion;
- implement and identify a nonnegative preconditioner approximation.

6.3 Workflow Gate B:

- sequence base pair -> FSC/cFAR -> NU evidence -> NU replay from one raw
  accumulation in both shared and distributed paths;
- remove the binary-envelope solvent prior and its automatic-mask bootstrap
  from the NU target path;
- reuse the compact evidence state for cutoff/local-resolution diagnostics and
  LP-set handoff rather than recomputing it after replay;
- leave the ordinary global-ML path bit-identical.

6.4 Science/cost Gates C and D:

- run only the two-way `P_tau` versus `Q_NU` ablation of §8;
- include uniform, heterogeneous, coarse-domain, weak-domain, and background-
  disorder fixtures plus real-data fixed-pose and refinement tests;
- judge at matched convergence, then measure the production-budget cost;
- do not begin Wilson work until Gate D is closed.

### Stage 7 — Wilson molecular precision — priority 2

Sequential only after Stage 6 Gate D. Start with a zero-mean shell-diagonal
Wilson spectrum behind one declared spectrum-source mechanism. Establish its
benefit against the accepted NU estimator as a separate experiment before any
combination is considered. Off-diagonal covariance, nonzero means, and local
Wilson/LocScale variants remain outside the first Wilson stage.

### Stage 8 — explicitly out of first scope

Singer/Gilles off-diagonal covariance; nonzero or phase-bearing prior means;
lagged LocScale surrogates; nonlinear/proximal priors; and any combined
NU/Wilson product-of-experts model. Each requires its own approved experiment
after the preceding estimator is independently understood.

## 11. The NU machinery as the prior infrastructure

The nonuniform-regularization machinery
(`doc/policies/nonuniform_filtering_policy.md`, `src/main/nu_filt/`) is the
mature evidence engine the prior now builds on directly. LocScale-2.0 derives
local confidence from a pseudoatomic reference; NU derives it from
cross-validated half-map prediction with an explicit noise model. The
2026-08-27 design joins the formerly separate envelope, local-resolution, and
graded-confidence ideas into one solver input.

### 11.1 Reuse without binary collapse

Reuse:

- the radially whitened symmetric cross-half Huber unary;
- the static candidate bank and accepted `nu_refine` extensions;
- spherical objective support and null-estimation discipline;
- ordered-label spatial regularization;
- selected-cutoff/local-resolution diagnostics and LP-set handoff.

Replace:

- coarsest-versus-best envelope margin as the sole prior statistic;
- binary MRF segmentation into molecular/solvent;
- component filtering, growth, cosine edge, and persistent envelope artifact;
- envelope complement and weighted-solvent-mean `Q_s`;
- a second post-replay NU analysis that can disagree with the prior state.

### 11.2 Required new evidence state

Add the null candidate and compact the full unary curve into ordered band
support/confidence before `dmats_mask` is released. The state is constructed
from current base halves between solve phases, not lagged from a previous
replay. It is shared read-only by the two half replays and by downstream NU
diagnostics. The base maps remain observable over the full spherical support,
so a weak domain excluded from the replay's useful bandwidth can still recover
evidence in a later refinement iteration.

Do not include the ML auxiliary replacement pair when constructing this state.
That pair is a post-hoc NU candidate under the current workflow and would let a
regularized result participate in estimating the precision that regularizes
it. `nu_refine` frontier evidence may participate only when the challenged
shell was evaluated from the unregularized base pair and accepted under the
recorded frontier contract.

### 11.3 Solver representation

The first operator is the compact multiscale frame of §5.3. Prefer broad bands
with an explicit normalization and a cheap positive preconditioner
approximation over an 8--16-label FFT stack. The selected cutoff is a
diagnostic/hand-off value; the solver consumes graded band confidence, not a
hard label mask. Measure the incremental value over post-hoc NU filtering:
the in-solve prior participates in deconvolution and Krylov conditioning,
whereas post-hoc filtering can only modify the finished estimate.

### 11.4 Feedback discipline

The envelope-specific lag-one shrinkage loop disappears, but resolution-field
feedback remains because later particle alignment consumes the NU replay map.
Track confidence/cutoff overlap and entropy across iterations, and retain the
omitted/weak/coarse-domain recovery tests. The ordinary base-pair FSC remains
the only resolution claim; shipped-pair correlation is diagnostic. Never feed
the NU evidence field into FSC solvent correction or replace the spherical NU
support with it.

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
