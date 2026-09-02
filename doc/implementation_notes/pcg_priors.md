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

**Outcome (2026-08-29): the Wilson question is answered and Stage 7 is
CLOSED.** `Q_W` was implemented, corrected twice by truth-judged failures,
validated, adjudicated against `Q_NU` on the 1WCM fixture, and REMOVED from
the code base the same day: `Q_NU` dominated at every shell even with `Q_W`
handed the ground-truth spectrum (S5.5 hypothesis confirmed). `Q_NU` is the
regularized estimator of this document; the Stage 7 sections below are the
complete experiment record. All forward-looking Wilson language elsewhere in
this document predates this outcome and is retained as plan history.

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
- **Solvent experiment REMOVED (2026-08-27)**: the binary-envelope solvent
  precision `Q_s`, its envelope loader, `pcg_solvent_lambda_rel`, the
  suppression/inflation convergence readout, and the harness strength-ladder
  sweep were removed from the codebase after the direct NU-evidence replay
  passed its first Gate C ablations (truth-judged win on 1WCM at both
  budgets; truth-free observables equal-or-better on bgal with the ML-replay
  stall fixed). §4 and the completed Stage 2--5 entries retain the algebra
  and measurements as the experiment record; `test=pcg_priors` was
  repurposed as the Q_NU Gate A suite. The abinitio3D automsk staging that
  the solvent bootstrap introduced is KEPT for its reference-masking
  registration benefit (Gate B streptavidin: base-pair 3.812 -> 3.191 A).
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
| Wilson molecular precision | CLOSED 2026-08-29 | adjudicated and removed; `Q_NU` dominates; §5.5 and the Stage 7 records |
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
decision superseded its binary molecular/solvent partition with the direct NU
precision of §5, and the code was REMOVED the same day once the NU
replacement passed its first Gate C ablations (see §1). The algebra, gates,
and real-data measurements below are retained so the experiment is not
repeated.

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

OUTCOME (2026-08-29): this experiment was run and its central question
answered in `Q_NU`'s favor -- see the Stage 7.1 records and the removal
decision. The hypothesis stated above (that the conservative NU estimator
already regularizes the unsupported degrees of freedom, leaving Wilson
extrapolation little to add) was CONFIRMED with `Q_W` at its best case.

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
  shared `regularize_state_half` and distributed
  `prepare_distributed_half_job`.
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
suppression coordinate. OPEN DECISION (superseded 2026-08-29, see the CLOSED
entry under Stage 6): lower PCG_SOLVENT_SUPP_OVER_PCT from
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

6.1 **IMPLEMENTED, build/runtime gate pending (2026-08-27).** Evidence contract,
no solver behavior change:

- `setup_nu_dmats(..., evidence_source='base_unfil')` opts into replay-evidence
  construction. The tagged path fingerprints the exact spherical-support values
  of both input halves, rejects auxiliary replacement candidates, and
  `build_nu_evidence_state` rechecks the fingerprint before compaction. Existing
  callers omit the tag and are unchanged; the current solvent prior remains in
  place until Stage 6.3 replaces it in the NU target workflow.
- The explicit null candidate is zero cross-half prediction: raw even is scored
  against zero odd prediction and raw odd against zero even prediction with the
  same radial whitening profile and Huber loss as the signal bank. It is
  smoothed at the discretized nominal-20-A scale, exactly matching its adjacent
  coarsest candidate. The raw zero loss has a systematic offset from smoothed
  predictors even under independent noise, and choosing the best signal label
  introduces a multiple-comparison advantage. The calibrated null therefore
  subtracts the robust median-plus-three-MAD offset of
  `C_zero-min(C_signal bank)` over the generous spherical support. This
  preserves the explicit cross-half null and lets truly coarse shared signal
  win whenever the coarse candidate's improvement exceeds the exact-bank null
  distribution. Null plus
  the full retained signal bank are regularized by a separate copy of the
  established ordered-label model, so the production NU filter map is not
  changed.
- The frozen `nu_evidence_state` has private packed storage and a copy-out API.
  It retains four nested coarse-to-fine support confidences (detail supported
  through 20/12/8/5 A), selected cutoff including cutoff zero for the null,
  normalized label entropy, spherical-support geometry, and scalar summaries.
  The geometry is sufficient to recreate the lexicographic packed-voxel order
  after mutable NU state is released. Confidence is a spatial-model
  softmax whose temperature is the median final best-versus-next energy gap for
  the exact bank. Provenance records the bank, ordering, whitening checksum,
  smoothing, temperature, and Potts scale; a content-derived FNV-1a identity is
  carried with the state. The null occupies one packed vector and reads the
  signal bank in place during compaction; it does not duplicate the full unary
  matrix. Validation enforces finiteness, `[0,1]` bounds, and monotone
  coarse-to-fine support.
- `simple_test_nu_envmask` now gates the null competitor on its independent-noise
  solvent/common-signal molecule: molecular coarse-band support must exceed
  solvent support, solvent null selection must exceed molecular null selection,
  all packed fields must satisfy the immutable-state contract, and the state
  must remain valid after mutable NU unary storage is released. Candidate-order
  mutation is guarded by a hard strict-order check; null calibration and exact
  provenance are emitted by `print_nu_evidence_summary`.
- Lightweight source validation completed: `git diff --check` and the Fortran
  source-index parser pass. The user-side build succeeded; runtime validation
  remains open before 6.1 is marked gate-complete.

First user-side execution reached compact-state construction after the NU bank
and ordered-label passes, then tripped the `[0,1]` immutable-state check. The
coarse support field sums all signal-label softmax probabilities and can exceed
one by single-precision accumulation roundoff. The derived support fields and
their summaries are now clipped at their mathematical probability bounds.
Softmax terms more than 80 calibrated energy units above the best state are set
to zero explicitly, avoiding the otherwise harmless IEEE underflow flag. Rerun
then showed that the uncalibrated zero predictor never won: molecular and
solvent null fractions were both zero, and coarse support saturated throughout
the support. This was a scientific null-model failure, not a test tolerance.
The first bias-calibration rerun still selected no null voxels: calibrating only
against the coarsest candidate did not account for another signal-bank member
winning solvent by chance. The calibration now uses the best exact-bank signal
cost, including that multiple-comparison advantage, without redefining the
coarsest signal label as solvent. Rerun of `simple_test_nu_envmask` is pending.

6.2 Operator Gate A:

- implement the normalized three- or four-band `Q_NU` factorization;
- prove adjoint identity, PSD, DC null, monotonicity, finite-difference
  gradient, and `P(H+Q_NU)P` composition;
- add the hard `P_tau`/`Q_NU` mutual-exclusion assertion;
- implement and identify a nonnegative preconditioner approximation.

6.2 **operator IMPLEMENTED, Gate A algebra tests pending (2026-08-27;
review-corrected the same day).** `set_nu_prior`/`apply_nu_precision` in
`simple_reconstructor_pcg`: `Q_NU = C (sum_b B_b^T W_b B_b) C` with
`B_b = crop o IFFT o M_b o FFT o pad`, `M_b` the disjoint radial 0/1 masks on
the padded Toeplitz lattice cut at the `NU_EVIDENCE_BAND_LIMITS` boundaries
(20/12/8/5 A; band 1 also absorbs everything coarser, band 4 everything
finer, padded DC and physically unaddressed points belong to no band), and
`C = I - 11^T/N` the native-box mean-centering projector applied on BOTH
sides. Review corrections to the first cut, recorded so the claims stay
honest:

- **Constant null mode.** Padded-DC exclusion alone is NOT sufficient: a
  native constant zero-pads to a box window with substantial non-DC padded
  content, so the first cut had `Q_NU(c 1) /= 0`. The explicit symmetric
  centering `C` on both sides now provides the EXACT constant null mode
  (Gate A must verify `Q_NU(c 1) = 0` to roundoff, not merely small).
- **Not a tight frame.** The pad/crop sandwich means the `B_b` are not
  orthogonal projectors and `sum_b B_b^T B_b < I` with cross-band leakage;
  refining or merging the band partition CHANGES the operator. The
  invariance claim is retracted. What survives, and is the declared
  normalization: each `B_b` is a contraction (restriction o projection o
  extension) and the `M_b` are disjoint, so
  `sum_b ||B_b x||^2 <= ||x||^2` and with `W in [0,1]`, `||Q_NU|| <= 1`.
  Gate A gains a partition-change test: uniform weights, two different band
  partitions, MEASURE the deviation and record it as the bank-refinement
  sensitivity rather than asserting zero.

`W_b = [p*(1-a_b)]^2` grades the lack-of-evidence weight by the soft support
exactly as the solvent weight did; each `B_b` is symmetric (restriction
adjoint to zero-extension, real even diagonal) and `C` is symmetric
idempotent, so `Q_NU` is symmetric PSD. Per-band results accumulate in real
space (two padded complex workspaces live). Strength:
`lambda_nu = pcg_nu_lambda_rel * s_data(D)`, derived in
`update_lambda_from_density` beside the ridge and solvent lambdas. Mutual
exclusion (R10) is enforced bidirectionally at set time in ALL pairs:
`set_ml_prior` rejects an attached `Q_NU` and vice versa, and `set_nu_prior`
and `set_solvent_prior` each reject the other (both attachment orders must be
mutated in Gate A). The declared nonnegative preconditioner approximation is
the support-mean band weight fused as a shell diagonal in
`finalize_density_accum`, mirroring the ML-prior fusion
(`nu_precond_shell_diag`). Deapodized-kernel-only attachment rides the
existing `assert_prior_attachment_mode`. Gate A algebra/mutation tests
(planned as a `test=pcg_priors` extension) are NOT yet written; the list now
includes the exact-constant-null check, the partition-change measurement, and
both mutual-exclusion attachment orders.

6.3 **workflow wiring IMPLEMENTED for both PCG paths, runtime gate pending
(2026-08-27).** `pcg_nu_lambda_rel` (default 0 = ordinary global-ML replay;
registered on `reconstruct3D` and `test=rec3D_backends`) selects the NU
replay: after both base solves and the FSC, `build_nu_replay_evidence` runs
`setup_nu_dmats(..., evidence_source='base_unfil')` ->
`optimize_nu_cutoff_finds` -> `build_nu_evidence_state` -> `cleanup_nu_filter`
on the current base half pair, prints the evidence block
(`pcg_replay_prior_mode=nu_evidence` + `pcg_nu_*` summary), and
`expand_nu_evidence_band_weights` (new, recreates the packed lexicographic
order from the frozen geometry alone; `w_b = 1 - a_b` inside the spherical
evidence support, 1 outside it) hands the solver its weights. Both half
replays share the one frozen state; `set_ml_prior`/`set_solvent_prior` are
never called in NU mode and no envelope is resolved, read, or written.
Review-added contracts (2026-08-27):

- **Readiness contract.** A valid compact state is necessary but not
  sufficient: `assert_nu_evidence_replay_ready` hard-errors before either
  replay when the explicit null wins less than `NU_EVIDENCE_MIN_NULL_FRAC`
  (provisional 1%, R9) of the generous spherical support -- the observed
  zero-null calibration failure must never attach silently -- or, after the
  first streptavidin run exposed the opposite failure, more than
  `NU_EVIDENCE_MAX_NULL_FRAC` (provisional 90%): a saturated null is equally
  a calibration failure, never a specimen property. The same assert
  now gates `simple_test_nu_envmask`, and that test must pass before the
  reconstruction harness results are trusted.
- **Explicit activation, no silent fallback.** `validate_nu_replay_request`
  (called from both PCG execution entries) hard-errors on a non-finite or
  negative strength and on a positive strength without the euclid ML replay
  (`ml_reg=yes`); `reconstruct3D` hard-errors on a positive strength with any
  non-pcg `rec_backend`. The harness deletes the key on its gridding leg.
- **Timing honesty.** The post-solve `get_nu_prior_stats` costs one full
  `Q_NU` application (~13 padded FFTs, material at small iteration budgets);
  it is timed separately, printed as `stats_overhead_s` on the
  `PCG NU REPLAY` line, and included in the shared path's replay total.

The
replay warm-starts from the previous-iteration shipped half when one exists
and otherwise cold-starts from the base solution (the closed-form shrinkage
init encodes the P_tau/Q_s optimum and is skipped). Per-half
`PCG NU REPLAY` lines report `lambda_eff` and `pcg_nu_prior_energy_final`;
the shipped-pair inflation crossings are measured in NU mode too. Distributed
parity: the same attach path runs in `prepare_distributed_half_job`; trailing
reconstruction and the trailing bootstrap are hard-errored with the NU replay
(the evidence contract requires the plain current-cohort base pair). The
solvent prior, its envelope bootstrap, and its default remain untouched
outside NU mode — removal from the target workflow is still owed once the NU
estimator passes its gates.

**Testable harness invocation (`test=rec3D_backends`).** A NU-replay run is a
single measurement (gates soft, as for any prior-active run), needs NO
envelope artifact, and forces `pcg_solvent_lambda_rel=0` explicitly (R10):

```text
simple_test_exec test=rec3D_backends projfile=<proj.simple> pgrp=<pgrp> \
    mskdiam=<D> pcg_nu_lambda_rel=<X> nthr=<N> [vol1=<truth.mrc> lp=<lp>]
```

against the strength-zero control (omit the key and use
`pcg_solvent_lambda_rel=0` for the pure `P_tau` reference, R5/Gate C run A)
and the solvent rungs already recorded. The execution directory carries a
`_nu<X>` token.

### Stage 6 first run: streptavidin NU replay at 0.1 (2026-08-27) — saturated-null calibration failure

First end-to-end harness run (10335 streptavidin, d2, box 128, mskdiam 80,
its=5, `pcg_nu_lambda_rel=0.1`). **Workflow invariants all held**: base pair
bit-identical between backends (3.518/3.119 A both legs, H2); one frozen
evidence identity for both half replays; solvent forced 0, no envelope
touched; `lambda_eff` 1.73e6 (vs the solvent prior's 1.53e6 on the same data
-- `s_data(D)` anchoring consistent); soft gates reported FAIL without
aborting, as designed. First H6 cost data: NU replay 26-31 s/half vs 6 s
base (~5x) plus 3-4 s `stats_overhead_s`.

**Scientific failure, opposite sign to the zero-null one:**
`pcg_nu_null_fraction=1.000000` with `uncertain_fraction=0` -- confidently
null everywhere including the molecule core, while the softmax still carried
28.6% coarse-signal mass (`band01=0.286`). `Q_NU` therefore degenerated into
a global detail penalty: shell amplitude ratio pcg/gridding declining
monotonically 0.9 (20 A) -> 0.33 (4 A), gated-band median 0.63 (FAIL at
<0.67). Diagnosis: the null offset `median + 3*MAD` of
`C_zero - min(C_signal)` over the whole support (median 0.37, MAD 0.18,
threshold 0.92) is a DETECTION threshold -- a signal label could only win
where its advantage exceeded the null population's 3-sigma -- and on a
molecule-dominated support the median itself is signal-contaminated. Right
shape for the retired binary envelope, wrong shape for a likelihood offset:
with overlapping components at 3.1 A it makes the null unbeatable. (It passed
the envmask fixture only because that molecule's gaps are hugely separated.)

**Calibration redesign (implemented same day):** the subtracted offset is now
the null component's CENTER only, estimated by the lower quartile of the gap
mixture (signal gaps sit strictly higher, so Q25 tracks the null component
even on a majority-molecule support); median/MAD stay recorded as
diagnostics, `null_offset=lower_quartile_center` in the provenance, and the
`+3*MAD` term is gone -- graded competition belongs to the softmax and
spatial consolidation to the ordered-label model. The readiness contract
gained the matching ceiling (`NU_EVIDENCE_MAX_NULL_FRAC=0.90`, provisional
R9): this run would now hard-error at evidence construction instead of
attaching a degenerate precision. **Envmask rerun under the center-only
offset: PASS (2026-08-27).** null_fraction 0.123 (inside the readiness
window); molecule coarse support 1.000 with zero null selections; solvent
support 0.829 with 18% null selections; band support monotone
0.884/0.884/0.804/0.291; envelope path untouched and passing. Watch item,
recorded not acted on: solvent coarse-band support 0.83 on an
independent-noise fixture means the Q25 offset is deliberately permissive
(~25% of gaps below it) -- weak coarse-band solvent penalty (W_1~0.17),
stronger fine-band penalty (W_4~0.7). Conservative direction for H3;
whether it regularizes ENOUGH is the S5.4 evidence-to-precision mapping
question the harness measures.

### Stage 6 first PASSING run: streptavidin NU replay at 0.1 (2026-08-27)

Identical invocation to the failed run; center-only offset (Q25 0.252 vs the
old detection threshold 0.92 on the same gap distribution: median 0.369, MAD
0.184). **PASS, all gates.** Evidence now graded and plausible:
null_fraction 0.189, band support 0.828/0.755/0.670/0.373 (monotone),
uncertain 3.1%. Observables:

- base pair pinned 3.518/3.119 (bit-identical to gridding and to the failed
  run -- H2/R5 hold);
- gated-band median amplitude ratio 1.066 (was 0.63 saturated-null), radial
  profile flat 1.00-1.09 inside 0.85x mask radius, centre-bin 0.79;
- shipped pair 3.267/3.049 -- 0.143 inflation vs base 2.2%, under the 5%
  bound, and equal to the ordinary P_tau shipped pair's 3.049 on this data;
- prior energy 7.1e9 (10x below the saturated run), ML RESID 1.7e-2 vs the
  P_tau-less base 4.0e-2 -- the precision also conditions the replay;
- known low-shell excess reappears mildly (k=1-2 ratio 1.75/1.50, inside the
  recorded historical 1.3-2.0 range);
- WATCH: near/beyond the band edge the NU replay carries more amplitude than
  the FSC-weighted gridding product (ratio rising 1.3 -> 7 over k=35-45,
  large beyond-band excess above k=46, outside the gated band). Expected
  structurally -- Q_NU has no per-shell Wiener rolloff, it suppresses only
  what local evidence disclaims (band-4 mean weight ~0.63) -- but whether
  that retained edge amplitude is signal or noise is precisely an R6
  question: only truth/held-out comparison decides it, never the half-map
  tables.

Next per Gate C: neutral-fixture run with ground truth (`vol1=` + `lp=`) and
the two-way P_tau-vs-Q_NU ablation at converged settings; then bgal (large
box, right-shifted dose-response) as the second dataset.

### Gate C first two-way ablation: 1WCM phantom, truth-judged — Q_NU wins (2026-08-27)

New user-built fixture (1WCM/RNA-Pol-II phantom, box 256, smpd 1.0, 2500
ptcls, c1, mskdiam 200, euclid+ml_reg, its=30 rtol=1e-4, lp=10 truth
comparison; an earlier B attempt at mskdiam=60 amputated the molecule inside
the support sphere -- negative low-k truth FSC, discarded, mask now matched).
Run A = `P_tau` control (`pcg_solvent_lambda_rel=0`), run B = `Q_NU` at 0.1.
Both PASS all gates; base pair bit-identical across A/B/gridding
(3.368/2.943 -- H2/R5). Evidence in B: null 0.126, bands
0.924/0.637/0.513/0.173, uncertain 0.

Map-to-truth FSC (the R6 decision variable; gridding column identical in
both runs):

| shell (A) | gridding | A: P_tau pcg | B: Q_NU pcg |
|---|---|---|---|
| 4.00 | 0.9427 | 0.9408 | 0.9543 |
| 3.20 | 0.7397 | 0.7317 | 0.7968 |
| 3.01 | 0.5931 | 0.5824 | 0.6748 |
| 2.91 | 0.5014 | 0.4951 | 0.5818 |
| 2.56 | 0.2079 | 0.1977 | 0.2844 |

Truth FSC=0.5 crossing: ~2.91 A (gridding), ~2.91 A (P_tau), ~2.84 A
(Q_NU). Q_NU matches everyone at low k (>=0.9988 through 5 A) and is
uniformly better than BOTH the P_tau replay and gridding from ~4 A outward
-- H1/H4/H5 supported on this fixture. Truth LS radial profile flat for
Q_NU (0.94-1.04 through bin 13) where gridding drifts to 1.3-1.5. Shipped
pair 3.241/2.844: 0.143 inflation vs base 3.4%, under the 5% bound.

Costs and watch items:

- **Convergence cost (H6).** The Q_NU replay ran 29-30 its at ~300 s/half
  (one half hit xtol at 29); the P_tau replay converged in 2 its / 6 s
  (warm shrinkage init ~= its optimum). ~50x converged-cost gap =
  Q_NU FFT stack x no closed-form initialization. Candidate mitigation for
  the production budget: a bandwise-shrinkage initial guess (scale each
  band by ~1/(1+lambda*W_b/rho) locally or shell-mean), the Q_NU analogue
  of `regularized_ml_initial_guess`. Production-budget (its=5) A/B still
  to be measured before any refinement integration.
- **Rim amplitude.** Q_NU map carries 1.5-2.0x gridding's |rho| in radial
  bins 10-13 (72-104 px) where P_tau's map declined; truth LS stays ~1
  and truth FSC is better, so it reads as retained real signal at the
  periphery rather than halo -- but verify against eroded/omitted-domain
  variants before trusting it generally.
- **Beyond-band retention.** amp_pcg/truth grows 20-80x for k>105 where
  truth FSC < 0.2: unshrunk beyond-band content survives at converged
  settings (S6 signature, outside the gated band). The evidence-mapped
  band-4 weight suppresses by lack of evidence, not by SSNR, so a
  Wiener-like rolloff is absent by design; postprocessing/FSC filtering
  owns the shipped-map rolloff as it does for the base maps.

### Gate C bgal ablation (no truth), its=30 rtol=1e-4 (2026-08-27)

Three-way on a fresh bgal registration (`bench.simple`, d2, box 256,
mskdiam 180, ~4.5k ptcls; base pair pins at 4.295/3.752 across ALL runs --
different registration than the recorded 4.132/3.709 ladder, internally
consistent): P_tau control, solvent Q_s at 0.1 (command-line slip, kept as a
free reference), Q_NU at 0.1. All three PASS all gates. Truth-free
observables (R6):

| observable | P_tau | Q_s 0.1 | Q_NU 0.1 |
|---|---|---|---|
| shipped 0.5/0.143 (A) | 4.185/3.667 | 4.132/3.667 | 4.132/3.627 |
| 0.143 inflation vs base | 2.3% | 2.3% | 3.3% (bound 5%) |
| gated-band ratio | 1.024 | 1.026 | 1.043 |
| agreement band | k=2-94 | k=2-94 | k=2-94 |
| radial norm min/max (in 0.85 mask) | 0.82/1.21 | 0.80/1.22 | **0.94/1.06** |
| centre-bin ratio | 0.76 | 0.76 | 0.85 |
| ML replay RESID / stop | 8.6e-2 xtol@12 | 1.4e-1 xtol@15 | **2.0e-2 maxits@30** |

Findings: (a) **the recorded bgal ML-replay stall is gone under Q_NU** --
the box-256 P_tau system stalls on step size at RESID ~0.09-0.14 while the
NU system reaches 2.0e-2, the second dataset where the shell-mean
preconditioner fusion conditions the NU replay BETTER than P_tau's; (b) the
Q_NU map has the flattest radial profile of the three and the best
centre-bin -- the rim erosion that dogged the solvent prior on bgal is
absent, consistent with the 1WCM rim-retention finding; (c) evidence sane
and mid-window (null 0.220, bands 0.763/0.664/0.576/0.169, uncertain ~0) --
the 0.90 ceiling was not approached despite the solvent-dominated support;
(d) shipped-pair inflation 3.3%, inside the bound. OPEN QUESTION (DEFERRED
2026-08-29, see the note at the end of this record) (the one
number that cannot be adjudicated without truth): the in-band amplitude
ratio pcg/gridding runs 1.5-2.0x through k=44-75 and higher toward the band
edge (in-band median 1.59 vs the control's 1.19) -- the expected signature
of no per-shell Wiener shrinkage on a low-SNR dataset (supported bands keep
base-level amplitude), and on 1WCM the same signature was truth-confirmed as
signal, but on bgal only an independent map-to-model comparison against the
deposited structure can decide it (R6). The postprocessing Wiener layer of
`nu_evidence_local_sharpening.md` is the designed consumer of exactly this
retained amplitude. Cost: 163 s/half at 30 its (production budget ~2 its
after warm start).

DEFERRED (2026-08-29): the map-to-model adjudication is parked by user
decision. It is not trivial in practice -- the obtained maps are not docked
to the deposited structure, so an honest comparison would require docking
and model building, which is not being undertaken now. The question stays
on record, not open-for-action. Confidence in reading the retained
amplitude as signal rests, for now, on the 1WCM truth-judged evidence
(where the identical signature was truth-confirmed, and where the 2026-08-29
Q_NU run was simultaneously amplitude-faithful and all-gates-passing);
this is supportive but is not a substitute for the real-data adjudication.
If the question is ever reopened, docking + model building against the
deposited bgal structure is the prerequisite; the LocScale-style Wiener
postprocessing layer remains parked behind it.

**Production-budget A/B (its=5, rtol=1e-3), same fixture (2026-08-27) --
H4 PASSES.** The Q_NU truth advantage survives the 5-iteration budget
essentially undiminished (truth FSC at 3.01 A: gridding 0.593, P_tau 0.580,
Q_NU 0.670 -- vs converged 0.593/0.582/0.675; at 2.91 A: 0.501/0.491/0.578;
at 2.44 A: 0.139/0.082/0.206). Both runs PASS all gates; base pair pinned
3.368/2.943 in all four runs; shipped-pair crossings identical to the
converged runs (Q_NU 3.241/2.844, inflation 3.4%). Notably the Q_NU replay
reaches RESID 8.4e-3 in 5 its (vs P_tau's 2.2e-2): the shell-mean
preconditioner fusion conditions the NU system well, and the converged-cost
gap was initialization distance, not ill-conditioning -- production budgets
do not suffer it. Cost at its=5: NU replay 56 s/half vs P_tau 13 s/half
(4.3x per replay; whole harness run 199 s vs 92 s including evidence build
and the 7 s stats overhead). Evidence stable across budgets (null 0.134 vs
0.126, bands 0.928/0.657/0.539/0.180). The bandwise-shrinkage init remains
worthwhile for converged/offline solves only.
Fallback if the quantile estimator proves fragile: calibrate the null on
noise-matched surrogates built from the half-map difference `(even-odd)`
(right noise spectrum, signal cancelled) -- costs a second bank pass and has
an anti-correlation subtlety (one independent difference realization only),
so it is recorded, not built.

6.3 Workflow Gate B:

- sequence base pair -> FSC/cFAR -> NU evidence -> NU replay from one raw
  accumulation in both shared and distributed paths;
- remove the binary-envelope solvent prior and its automatic-mask bootstrap
  from the NU target path;
- reuse the compact evidence state for cutoff/local-resolution diagnostics and
  LP-set handoff rather than recomputing it after replay;
- leave the ordinary global-ML path bit-identical.

6.3 **refinement integration + solvent removal DONE (2026-08-27), runtime
gate pending.** `pcg_nu_lambda_rel` is registered on `refine3D` and
`abinitio3D` (activation-gated on `rec_backend=pcg`) and forwarded through
`apply_refine3D_reconstruction_controls`, exactly the `maxits_pcg`/`rtol`
route the solvent key used -- an abinitio3D run opts in with
`rec_backend=pcg pcg_nu_lambda_rel=0.1` and every stage's regularized
replays run the NU precision with the established cross-iteration warm
start. The solvent prior is fully removed: `Q_s` and its stats from
`simple_reconstructor_pcg`; the envelope loader, suppression reference,
stats reporter/persistence, and the envelope-flattening part of the
regularized init from the strategy (both paths); the suppression/inflation
readout, guidance thresholds, and stats reader from `simple_convergence`;
`pcg_solvent_lambda_rel` from parameters/parse and all four UIs;
`PCG_SOLVENT_STATS_FILE` from the fname defs; the strength-ladder sweep,
envelope requirement, envelope symlinking, and solvent summary fields from
`test=rec3D_backends` (now always a single comparison; NU runs get soft
gates). `test=pcg_priors` is repurposed as the Gate A Q_NU suite: adjoint,
PSD, EXACT constant null (the centering mutation target),
zero-action-under-full-support, monotonicity under evidence withdrawal, FD
gradient, the MEASURED 4-band-vs-2-band partition sensitivity, masked
composition, and priored solve parity across the shared-vs-nparts reduction
seam. The abinitio3D automsk-from-first-NU-stage staging is kept for the
reference-masking registration benefit; the shipped-pair crossing
measurement is kept as the NU over-regularization diagnostic (in-strategy
print only; the solvent stats file and convergence guidance are gone).
Build + `test=pcg_priors` + `test=rec3D_backends` reruns pending user-side;
first abinitio3D run with the NU prior is the outstanding refinement gate.

### First NU-replay abinitio3D: bgal (2026-08-28) — refinement gate observations

`rec_backend=pcg pcg_nu_lambda_rel=0.1`, distributed, stages 4-6+ observed
(box crop 130 -> 140, lp 10.2 -> 9.3 -> 4.66 promoted by the NU filter):

- **End-to-end NU refinement loop works.** Every iteration: base solves (2
  its, ~1 s/half) -> fresh evidence from the current base pair (distinct
  identity per iteration, readiness passing) -> warm start from the previous
  shipped half -> 2-it Q_NU replay (~4.5-7 s/half + ~1.1-1.8 s stats
  overhead). Production-viable cost: reconstruction remains a small fraction
  of the alignment wall time.
- **Evidence trajectory is stable and physically sensible** (the S6.3
  field-overlap tracking item, answered by inspection): null pinned at
  0.229-0.232 across every iteration and stage; coarse bands constant
  (band01 ~0.764); band04 support 0.000 through stage 5 (crop Nyquist above
  5 A -- correctly no fine evidence), then 0.086 -> 0.146 -> 0.161 as stage
  6 extends the band and resolution genuinely develops. No oscillation, no
  runaway. Warm-started ML RESID 0.056-0.075 at 2 its (base 0.13).
- **Resolution progression**: 5.44 -> 5.02 -> 4.66 A at 0.143 through
  stages 4-6, cFAR 0.64 -> 0.75.
- **The solvent-era stage-6 overlap collapse is largely CURED**: at NU/lp
  engagement (lp 9.3 -> 4.66 + switch to envelope-masked references) the
  orientation overlap dips 0.92 -> 0.80 -> 0.60 and recovers to 0.86 within
  two iterations -- versus the recorded 0.70 -> 0.04 collapse with the
  solvent prior. The NU replay map is evidently a much better matching
  reference at promotion.
- **Kept automsk staging coexists correctly**: stage 6 first iteration
  generates the envelope (occupancy 29.6%, one component), later iterations
  consume it for reference masking only; no solvent-prior lines anywhere.
- **Remaining S6.3 sub-item now visible as measured cost**: the NU analysis
  runs TWICE per stage-6 iteration -- once for the replay evidence, once for
  the volassemble NU filter/envelope postprocess (identical whitening and
  ordered-label logs back to back). Sharing the compact evidence state with
  the postprocess remains the outstanding dedup, now with a concrete
  per-iteration price attached.

### DECISION (2026-08-28): Q_NU default-on in NU mode; post-hoc NU filter and pcg automasking retired from the NU-replay path

Following the first NU-replay abinitio3D (above) and the user-side
observation that the Q_NU shipped map and the post-hoc NU-filtered map are
nearly indistinguishable, three policy changes landed together:

1. **Default-on.** `pcg_nu_lambda_rel` dynamically defaults to 0.1 whenever
   `rec_backend=pcg`, the NU machinery is active (`l_nonuniform`), and the
   euclid ML replay runs -- i.e. the NU-filtered stages of abinitio3D/refine3D
   run the Q_NU replay with no flags. An explicit 0 restores the P_tau
   replay (the R10 A/B control); the plain-`reconstruct3D`/harness default
   (no `filt_mode`) stays 0, preserving R5 for the recorded baselines.
2. **No post-hoc NU filtering with Q_NU in the solve.**
   `filter_pcg_nonuniform_maps` returns early when the NU replay is active:
   no second NU analysis (S6.3 dedup closed), no `_nu_filt`/`_nu_locres`
   products, no evidence envelope (envelope clause amended 2026-08-29, see
   the automsk decision below). The LP-set matching handoff survives,
   now derived from the frozen replay evidence itself (finest evidenced
   local cutoff per state, threaded
   `build_nu_replay_evidence -> execute_rec3D_pcg_distributed_master ->
   filter_pcg_nonuniform_maps`, same `set_all2single('lp',...)` contract).
   Matchers fall back to the regular (Q_NU-regularized) references, which
   is the point. The full postprocess remains for the P_tau/gridding paths.
3. **No forced automasking in the pcg path of abinitio3D.** The staging
   existed to bootstrap the solvent envelope and switch to envelope-masked
   references; with the prior gone and the Q_NU map already locally
   regularized, automasking follows the same explicit user control as
   everywhere else. (The Gate B registration-benefit observation attributed
   to envelope-masked references is superseded by the user's read that the
   NU-replay reference makes it unnecessary; if registration regresses on a
   future run, this is the first knob to revisit -- and since 2026-08-29 that
   knob works: explicit `automsk=yes` regenerates the envelope from the
   frozen replay evidence without resurrecting the post-hoc filter, see the
   automsk decision below.)

Verification pending: an abinitio3D rerun with zero prior flags showing the
dynamic default engaging at the first NU stage, single NU analysis per
iteration, evidence-derived lp promotion, and no envelope/automask lines.
CLOSED 2026-08-28/29 by the three-dataset run record below.

### Stage 6 record: zero-prior-flag abinitio3D verification runs (2026-08-28/29, msp1 + embb + exp_gate)

Three abinitio3D runs with `rec_backend=pcg` and NO prior flags
(commit 6d06a1fa; logs under `Processing/pcg_integration/`), closing the
default-on decision's pending verification. Observed on the completed runs
(msp1, 5947 s; embb, 4228 s; exp_gate consistent while still in flight):

- dynamic default engaged at every NU stage (`Q_NU on, LAMBDA_REL 0.1`),
  early non-NU stages correctly `off`; evidence built once per state per
  volassemble from the FSC half pair (`source=base_unfil`), no second NU
  analysis, no `_nu_filt`/`_nu_locres`/envelope/automask lines;
- evidence-derived matching low-pass promoted every cycle
  (`set_all2single('lp',...)` values tracking the base FSC=0.5 region);
- no over-regularization: shipped-pair FSC=0.143 crossings track the base
  crossings essentially exactly (msp1 3.91-3.96 vs 3.86-3.96 A; embb both
  4.00 A), never pulling finer;
- every solve stopped on `fixed_iterations` (2 its) with residuals of a few
  percent and near-identical even/odd values; halfset cohorts balanced;
  trailing `F == U` with the expected stage fraction ramps; ML warm starts
  firing throughout;
- stable, slowly rising evidence support (band01 ~0.77 -> 0.78); final
  resolutions Nyquist-limited at the stage crops (msp1 3.86 A, embb 4.00 A);
  refinement cFAR ~0.80 (msp1) with a 0.51 final full-crop bootstrap value
  (real directional anisotropy, not a mechanism issue); B-factors -103/-58.
- convergence contrast worth keeping an eye on, not a gate item: embb
  converged classically (orientation overlap 0.35 -> 0.87, ~1.6 deg mean
  angular distance) while msp1/exp_gate plateaued at overlap ~0.10-0.22
  with within-stage climbs and no collapse signature (the periodic dips
  align with fill-in sampling iterations).

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

DECISION (2026-08-29): Wilson development STARTS now, ahead of formal Gate D
closure, by user direction -- the Stage 6 verification runs and the Gate C
ablations already on record justify overlapping the engineering. The
sequencing constraint that matters is preserved unchanged: Wilson is
compared against the NU estimator as a separate experiment, is never
combined with it (or with P_tau) in this stage, is off by default behind an
explicit opt-in, and cannot be considered for any default until Stage 6
Gate D and its own gates close. See the Stage 7 implementation plan below.

### Stage 7.1 — Q_W operator and workflow wiring (historical record: IMPLEMENTED, VALIDATED, ADJUDICATED, and REMOVED 2026-08-29 — Q_NU dominates, see the Gate C record and the removal decision)

**Form.** `Q_W = F* diag(q_W(k)) F`, zero-mean and shell-diagonal, exactly
the first form S5.5 prescribes (anchoring corrected twice by the first two
runs, see the S7.1 run records below):

```text
q_W(sh) = lambda_rel * Dbar(sh) * min(WILSON_PREC_CAP, s(k_hp) / s(sh))
WILSON_PREC_CAP = 1e2
```

`s` is the Wilson spectrum SHAPE from the declared source, anchored at the
first prior-active shell `k_hp`; `Dbar(sh)` is the shell-mean raw data-only
density -- the same per-shell statistic `P_tau` divides by `tau*SSNR`. The
per-shell D anchoring makes the prior/data ratio exactly
`lambda_rel * min(CAP, s_ref/s(sh))` at EVERY shell, bounded in
`[lambda_rel, CAP*lambda_rel]` and independent of the CTF/sampling decay
of D -- the bounded-dynamic-range discipline `Q_NU` established (bounded
operator, one strength knob), expressed in `P_tau`'s per-shell data
convention. The spectrum supplies shape only. Shells coarser than the `hp`
limit get no prior, mirroring the `P_tau` low-frequency no-prior guard.

**No parallel Fourier array.** The builder
(`build_wilson_prior_from_spectrum`) fills the SAME calibrated `ml_prior`
diagonal `P_tau` uses, with the identical padded-radius -> native-shell
mapping. Application in the replay operator, kernel-diagonal fusion, the
preconditioner contribution, and the `get_ml_prior_stats` readout are
therefore shared verbatim -- `Q_W` is a different *builder* behind the one
diagonal, and three-way mode exclusion (`P_tau` | `Q_NU` | `Q_W`) is both
structural (one diagonal) and asserted in all three setters
(`set_ml_prior`/`set_nu_prior`/`set_wilson_prior`, R10 convention).

**Declared spectrum source (7.1).** An explicit reference volume,
`pcg_wilson_vol`: read via `read_and_crop` to the solve grid, rotationally
averaged power spectrum (`image%spectrum('power')`), loaded once per
execution and shared by both halves like the `P_tau` FSC prior. For the
phantom gates this is the ground-truth map -- the best-case Wilson prior,
so the first experiment answers "does a CORRECT molecular spectrum help"
before any estimation machinery is built. Estimated sources (Wilson-line
fit of the supported band with extrapolation) and analytic
composition-based spectra are later sub-stages (7.2+), each behind the same
single declared-source mechanism.

**Activation contract.** Explicit opt-in only: `pcg_wilson_lambda_rel > 0`
(no dynamic default) requires `rec_backend=pcg`, the euclid ML replay,
`pcg_wilson_vol`, and ordinary (non-NU) mode -- `Q_W` + NU filtering
hard-errors, which keeps the first-stage comparison clean (no post-hoc NU
filtering of Wilson maps, no dynamic `Q_NU` default interference) and makes
`Q_NU`/`Q_W` mutual exclusion follow from the mode split; both exclusions
are ALSO hard-asserted at parameter validation and in the reconstructor.
A positive strength on any other backend hard-errors (explicit-activation
contract, S6.2). Both halves cold-start from the base solution (the
shrinkage initial guess encodes the `P_tau` optimum; no Wilson closed form
yet); the own-half warm start applies as everywhere.

**Wiring inventory.** `set_wilson_prior`/`build_wilson_prior_from_spectrum`
+ state in `simple_reconstructor_pcg`; `load_wilson_spectrum` + replay
branches in both PCG execution paths of `simple_rec3D_pcg_strategy`;
`pcg_wilson_lambda_rel`/`pcg_wilson_vol` in parameters (declaration, parse,
validation phase); UI registration for `reconstruct3D`, `refine3D`, and the
`rec3D_backends` test harness; backend guard in the `reconstruct3D`
commander. The `rec3D_backends` harness handles the Wilson keys like the
NU key: validated up front (ml_reg required, Q_NU exclusion), soft gates
when the prior is active (it legitimately moves the pcg leg away from
gridding), `_wil<strength>` execution-directory tag, `pcg_wilson_vol`
made absolute before the chdir, and both keys deleted from the gridding
leg.

**Stage 7.1 gates (as planned; overtaken by events -- Gate C was run first
and concluded the stage, so the Gate A/B formalization below was never
executed and is moot after the removal).**

- Gate A algebra: `Q_W` diagonal positivity and shell-constancy on the
  padded lattice; the three-way exclusion mutation tests (each setter pair
  rejected in both attachment orders); prior-stats readout sanity
  (`prior_to_khat` ratios finite, hp shells zero).
- Gate B invariants: artifact set identical to the `P_tau` replay
  (`Q_W` changes only the diagonal); FSC/resolution reporting unchanged;
  the spectrum-source log line present exactly once per execution.
- Gate C: truth-judged 1WCM two-way ablations `Q_W` vs `P_tau` and `Q_W`
  vs `Q_NU` at matched convergence, with the truth spectrum as source
  (best case); then bgal with the deposited-model spectrum. Judge on the
  same criteria as the Stage 6 ablations; the S5.5 hypothesis to test is
  that the conservative NU estimator already captures most of the benefit
  and Wilson extrapolation adds little on data the experiment observed.

### Stage 7.1 first run: 1WCM truth-judged (2026-08-29) — spectrum-peak anchoring failure

First `rec3D_backends` run with the truth map as spectrum source
(`pcg_wilson_lambda_rel=0.1 pcg_wilson_vol=1WCM.mrc`, euclid+mlreg, its=5,
box 256 smpd 1.0). The harness worked as intended: both legs ran, soft
gates reported, and the truth tables adjudicated. Result: FAILURE of the
initial `Q_W` anchoring, cleanly diagnosed.

- The base pcg solves were healthy (rel residual ~1.2e-2, base-pair
  FSC=0.143 at 2.94 vs gridding 2.91 A, cFAR 0.89 vs 0.90) -- the failure
  was confined to the `Q_W` replay, as the mode split predicts.
- The `Q_W` replay solves left rel residual ~42-45 after 5 iterations
  (vs ~3e-2 for `Q_NU`/`P_tau` replays): the system was unsolvably stiff.
  The shipped map was a transient-dominated CG iterate: fsc(truth,pcg)
  NEGATIVE through the mid band (to -0.74), amplitude collapsed to a flat
  floor mid-band and rising junk toward Nyquist, shipped-pair FSC
  crossings far coarser than base.
- Cause: the initial form max-normalized the spectrum shape and floored it
  at 1e-4. A real map's power spectrum is dominated by the lowest shells
  by orders of magnitude, so ESSENTIALLY EVERY prior-active shell sat at
  the floor and received precision ~1e3-1e4 x data_scale: the prior
  crushed the data term band-wide.
- Fix (implemented): anchor the shape at the first prior-active shell
  (`k_hp`) and cap the precision dynamic range at `WILSON_PREC_CAP = 1e2`
  -- `q_W in [lambda_w, 1e2*lambda_w]` over the active band, preserving
  the Wilson decay shape where it is resolved and saturating beyond. This
  is the same bounded-operator discipline as `Q_NU`. Rerun below.

### Stage 7.1 second run: 1WCM truth-judged (2026-08-29) — anchoring fix verified; global-data_scale over-shrinkage identified

Identical invocation after the hp-anchoring fix. Solvability restored and
the first genuine `Q_W` measurement obtained:

- `Q_W` replay residuals 1.2e-1 (from 42-45), truth-FSC positive across
  the band, shipped map a real regularized estimate.
- LOW-BAND WIN against the ML-regularized gridding leg:
  fsc(truth,pcg) >= fsc(truth,gridding) through k=2-10, decisively at
  k=3-5 (0.998 vs 0.956 at k=3) -- the molecular spectrum is informative
  exactly where Wilson statistics say it should be.
- Progressive high-band over-shrinkage: pcg truth-FSC falls behind
  gridding from mid-band (k=80: 0.59 vs 0.74), amp_pcg/truth decays to
  ~0.29 by k~75, and the shipped-pair crossing moved COARSER than base
  (3.08 vs 2.94 A) -- the over-regularization direction.
- Diagnosis: anchoring `lambda_w` to the GLOBAL data_scale while the
  data term `Dbar(sh)` decays with k (CTF + sampling) inflates the
  effective prior/data ratio k-dependently by `data_scale/Dbar(sh)` on
  top of the spectrum ratio -- a mis-anchoring no single lambda can
  compensate.
- Fix (implemented): per-shell D anchoring,
  `q_W(sh) = lambda_rel * Dbar(sh) * min(CAP, s_ref/s(sh))`, mirroring
  `P_tau`'s convention exactly; the prior/data ratio is now bounded in
  `[lambda_rel, CAP*lambda_rel]` at every shell. Rerun below.

### Stage 7.1 third run: 1WCM truth-judged (2026-08-29) — per-shell anchoring VALIDATED; Q_W wins both band ends

Identical invocation after the per-shell D anchoring. The operator form is
now measured-correct and the first Wilson science result is on record:

- Solve health at `Q_NU` class: replay residuals 3.0e-2 (from 1.2e-1);
  shipped-pair FSC crossings IDENTICAL to base (3.368/2.909 A) -- zero
  over-regularization signature at lambda=0.1.
- LOW-BAND WIN retained: fsc(truth,pcg) >= gridding through k=2-9
  (0.9993 vs 0.9908 at k=2; 0.9981 vs 0.9556 at k=3).
- HIGH-BAND WIN gained: from k~79 to the band edge pcg matches or beats
  the ML-regularized gridding reference against truth (k=80: 0.743 vs
  0.740; k=90: 0.455 vs 0.440; k=100: 0.225 vs 0.209; k=110: 0.100 vs
  0.058) -- the molecular prior paying off exactly where Wilson statistics
  extrapolate into weakly measured shells. Run 2 lost this entire region.
- Remaining: a small systematic mid-band truth-FSC deficit (k~30-75,
  0.005-0.02, e.g. 0.964 vs 0.974 at k=50) with the amplitude ratio
  plateauing ~0.86 -- honest lambda=0.1 shrinkage, the lambda sweep's
  target. High-k amplitude rise beyond the resolved band (rim shells)
  persists but with better truth-FSC than gridding there.

Next: the lambda sweep (0.01/0.03/0.1) to see whether a weaker prior
closes the mid-band deficit while keeping both band-end wins, and the
`Q_NU` two-way on the same project (`pcg_nu_lambda_rel=0.1` in place of
the Wilson keys) -- the first direct `Q_W` vs `Q_NU` adjudication. These
are the Gate C entry points. Both executed below.

### Gate C: Q_W lambda sweep + first Q_W vs Q_NU two-way, 1WCM truth-judged (2026-08-29) — Q_NU DOMINATES; S5.5 hypothesis CONFIRMED

Same fixture and budget throughout (euclid+mlreg, its=5, truth spectrum as
the Q_W source -- Wilson's best case).

**Q_W lambda sweep (0.1 / 0.03 / 0.01):**

- The mid-band truth-FSC deficit against the gridding reference shrinks
  monotonically with lambda (k=50: 0.964 / 0.967 / 0.972 vs gridding
  0.974) while the high-band win is essentially lambda-INSENSITIVE
  (k=110: ~0.10 at all three strengths vs gridding 0.058), and the
  low-band win is common to all pcg runs. Within the family,
  lambda~0.01 dominates: near-zero mid-band cost, full high-band gain,
  replay residuals 1.0e-2. The high-band gain evidently comes from
  suppressing noise in weakly measured shells at all, not from the
  precise Wilson strength -- a first hint that the spectrum SHAPE is not
  the operative ingredient.
- Amplitude: mid-band amp/truth ~0.86 / ~2.3 / ~4.6 (vs gridding ~9) --
  the single global knob trades amplitude fidelity against noise
  suppression band-wide, as a shell-diagonal must.

**Q_W vs Q_NU two-way (pcg_nu_lambda_rel=0.1):**

`Q_NU` beats the gridding reference AND every `Q_W` run at EVERY shell,
including the high band that was Wilson's presumptive niche (truth-FSC
k=50: 0.978; k=70: 0.925; k=80: 0.793; k=90: 0.523; k=100: 0.284; k=110:
0.137 -- vs gridding 0.974/0.905/0.740/0.440/0.209/0.058 and best-Q_W
0.972/0.901/0.738/0.445/0.219/0.099). It does so while remaining
amplitude-faithful: amp/truth ~9 band-wide (no shrinkage), the flattest
radial LS profile of any run including gridding (0.89-1.02 vs gridding's
1.14-1.22 over-amplification), and the harness passed ALL gates --
including the gridding-parity amplitude gates every other prior run
soft-failed. Shipped-pair inflation modest (3.24/2.84 vs base 3.37/2.94),
suppression ~36%, replay residual 8.4e-3. Cost is the one `Q_W`
advantage: ~60 s/half (15 s evidence overhead) vs ~13 s.

**Verdict.** The S5.5 hypothesis is CONFIRMED on this fixture: the
conservative NU estimator captures more than the full Wilson benefit --
even with `Q_W` handed the ground-truth spectrum -- because the evidence
field regularizes locally and anisotropically exactly the degrees of
freedom the half-map evidence does not support, which no shell-diagonal
can express. `Q_W` is NOT a candidate to displace `Q_NU` and no further
solo `Q_W` development is warranted on this evidence (removal decision
below). Remaining bookkeeping if ever needed: a bare pcg run (no prior
keys) on this fixture for the clean `Q_W` vs `P_tau` two-way; not
required for the verdict above.

### DECISION (2026-08-29): Q_W REMOVED from the code base

Following the Gate C verdict, the Wilson estimator was removed the same
day it was implemented, by user direction and in keeping with this
document's discipline: retired estimators live on as experiment records,
not as code. Removed: `set_wilson_prior`/`build_wilson_prior_from_spectrum`
and all Wilson state/constants from `simple_reconstructor_pcg`; the
spectrum loader and replay branches from both PCG execution paths;
`pcg_wilson_lambda_rel`/`pcg_wilson_vol` from parameters, parsing,
validation, and the `reconstruct3D`/`refine3D`/`rec3D_backends` UI; the
backend guard from the `reconstruct3D` commander; and the Wilson key
handling from the `rec3D_backends` harness. The three-way mode assertion
reverts to the two-way `P_tau`/`Q_NU` exclusion. The S7.1 records above
(operator form, both anchoring corrections, the sweep, and the two-way)
are the complete archaeological record; if a Stage 8 combination
experiment is ever approved, the validated per-shell-D-anchored form
documented here is the reference implementation to resurrect. Wilson
remains formally priority 2 in name only; in practice Stage 7 is CLOSED.

### Stage 8 — explicitly out of first scope

Singer/Gilles off-diagonal covariance; nonzero or phase-bearing prior means;
lagged LocScale surrogates; nonlinear/proximal priors; and any combined
NU/Wilson product-of-experts model. Each requires its own approved experiment
after the preceding estimator is independently understood.

Related but separate (POSTPROCESSING, not a solver prior, so the in-solve
LocScale risks above do not apply): model-free LocScale-style local
sharpening driven by the same frozen NU evidence state is specified in
`nu_evidence_local_sharpening.md` — parked until this document's Gate C/D
program completes.

### DECISION (2026-08-28): NU-replay firing readout in the convergence report

With `Q_NU` default-on in NU mode, refinement needs the same in-run answer the
retired solvent prior had: "is the prior firing, and is
`pcg_nu_lambda_rel` right?" — with no harness and no ground truth. The readout
is the **NU prior-energy suppression percentage**, per half:

    supp% = 100 * (1 - sqrt(E_NU(x_ML) / E_NU(x_base)))

where `E_NU(x) = (lambda_nu/2) x^T Q_NU x` (`get_nu_prior_stats`; the
`lambda_nu` factor cancels in the ratio) and `x_base` is the unregularized
base solution of the SAME half — the replay's own reference: with `P_tau`
absent in NU mode and the replay cold-started from the base solution, a
vanishing `pcg_nu_lambda_rel` reproduces `x_base`, so an inert prior reads
~0% with no shrinkage referencing needed (unlike the solvent readout, whose
reference had to factor out the Wiener component of the `P_tau` replay). The
amplitude-domain square root mirrors the solvent readout's rms ratio. Cost:
one extra full `Q_NU` application per half (~13 padded FFTs), timed into the
existing `stats_overhead_s` diagnostic.

Plumbing mirrors the retired solvent readout: `report_nu_solve_stats` prints
`pcg_nu_prior_energy_final/base` and `pcg_nu_suppression_pct` per half; both
PCG execution paths average even/odd per state and persist
`simple_pcg_nu_stats.txt` (delete-then-rewrite each volassemble, so a
skipped replay leaves no stale values), including the shipped-pair FSC=0.143
crossing per state (the over-regularization diagnostic, never a resolution
claim); `check_conv3D` prints `% PRIOR ENERGY SUPPRESSED (PCG NU REPLAY)` and
`SHIPPED-PAIR FSC=0.143` with the other iteration stats, writes them to
`simple_stats.txt` (`PCG_NU_SUPPRESSION_PCT`, `PCG_NU_SHIP0143`,
`PCG_NU_LAMBDA_REL` + per-state keys), and advises:

- **< 5%**: prior inert — increase `pcg_nu_lambda_rel` (~3x)
- **5-60%**: nominal — keep the strength
- **> 60%**: over-regularization risk — decrease (~3x)

Thresholds are provisional constants in `simple_convergence`
(`PCG_NU_SUPP_INERT/OVER_PCT`), inherited from the solvent readout's bounds;
recalibrate with a `pcg_nu_lambda_rel` strength ladder now that every run
prints the readout, reading the suppression % against the shipped-pair
crossing and the evidence trajectory.

### DECISION (2026-08-28): NU replay supports trailing reconstruction (evidence pair = FSC pair)

The hard block on `Q_NU` + `trail_rec` is lifted; real-data refinement needs
both. The governing contract is **the evidence pair is always the FSC pair of
the volassemble**:

- **Trailing, chain present**: the base solves are the full-mass blended
  chain solutions — the very statistics the ML replay re-reads
  (`add_raw_accum_weighted` of the same chain) — so current-iteration
  evidence from that pair (`source=base_unfil`) satisfies the "evidence from
  the pair the replay reuses" requirement. The original block's concern (the
  base pair is not the plain current-cohort pair) is resolved by matching the
  evidence to the replayed statistics, not to the cohort.
- **Trailing bootstrap (no chain yet)**: the FSC already comes lag-one from
  the previous iteration's shipped half pair; the evidence comes from that
  same pair (`source=previous_shipped`, a new allowed evidence source next
  to `base_unfil`). Lag-one evidence is precedented by the retired solvent
  envelope's one-iteration lag and is exactly as trustworthy as the
  bootstrap FSC that `P_tau` would be built from. The current-cohort
  fractional pair is NOT used: at small update fractions it can legitimately
  fail null calibration and would hard-stop real runs.

The accumulator arithmetic (chain blend weights, `scale_raw_accum`,
fixed-order reduction) is the `test=pcg_frac_update`-gated path
(`validate_rec3D_pcg_fractional_updates`) and is untouched: `Q_NU` attaches
after accumulation and before `end_accum`, exactly where `set_ml_prior`
does, so the prior is orthogonal to the tested trailing contracts. The
firing readout's reference remains the base solution of the same half (in
trailing, the blended base solve the replay warm-starts from), and the
LP-set handoff derives from whichever evidence pair was used. Shared-memory
trailing remains unsupported for PCG generally (pre-existing
accumulator-domain restriction, unrelated to the prior).

### DECISION (2026-08-29): automsk=yes regenerates the evidence envelope under the Q_NU replay

Amends item 2 of the default-on decision: with `automsk` enabled, the NU
evidence envelope IS produced on the replay path -- regenerated inside
`build_nu_replay_evidence`, between `optimize_nu_cutoff_finds` and
`cleanup_nu_filter`, while the raw per-voxel evidence margins are live. The
matching-reference envelope therefore derives from the same frozen evidence
as the `Q_NU` precision: one analysis, two consumers, and still no second NU
pass. Cadence and artifact naming follow the same `plan_state_postprocess`
contract as the post-hoc paths (missing/incompatible mask, `startit`,
`AMSK_FREQ`). All three regeneration sites (gridding volassemble, PCG
post-hoc, PCG replay) now share the single producer
`write_nu_evidence_envmask` in `simple_nu_filter`; the replay envelope is
built from the static candidate bank (the replay analysis runs no
high-resolution shell extension), a deliberate, slightly shallower evidence
basis than a gridding-path envelope. This also arms the default-on
decision's item-3 fallback: if reference registration regresses without
envelope masking, explicit `automsk=yes` is now a functional knob on the
replay path.

### DECISION (2026-08-29): refine3D_auto joins the pcg bypass; NU shell extension declared obsolete under Q_NU

`refine3D_auto` was the last workflow defaulting into the competition NU
machinery on the pcg backend: its unconditional `nu_refine=yes` default
hard-errored at the first volassemble, and its init-volume bootstrap ran the
competition prefilter (including the shell challenger) on the startup
halves. Both are fixed, mirroring the abinitio3D stage policy: `nu_refine`
defaults to `no` when `rec_backend=pcg` (explicit `nu_refine=yes` + pcg
remains a hard error), and the startup NU prefilter is bypassed under pcg --
raw-pair validation and its reconstruct-startup fallback stay
backend-neutral, but no `_nu_filt` startup references are produced and the
first matcher pass uses the raw native E/O references, exactly as every
later Q_NU iteration feeds the matcher. Doctrine going forward:
`rec_backend=pcg` with the `Q_NU` prior makes NU high-resolution shell
extension obsolete; `nu_refine` survives only on the gridding/competition
path pending its retirement (next section).

### Gate C/D comparability note (2026-08-28 mask unification)

The gridding backend's spherical FSC mask was unified to `params%msk_crop`
(previously the broad rim radius), so gridding and PCG FSC/cFAR are now
directly comparable -- but gridding FSC-derived numbers recorded BEFORE the
unification (including baselines quoted in the Stage 5/6 records above) are
not comparable to post-unification gridding runs. Remaining Gate C/D
ablations must use same-policy baselines: rerun the gridding arm at the
current mask policy, or compare against the PCG backend directly.

### CLOSED (2026-08-29): solvent-era suppression-threshold OPEN DECISION superseded

The open decision under the Stage 5 streptavidin record (lower
`PCG_SOLVENT_SUPP_OVER_PCT` from 60 to ~45-50; move the solvent default
0.1 -> 3e-2) is superseded: the solvent prior is scientifically retired and
`Q_NU` is the default regularized estimator, so the solvent suppression
guidance thresholds will leave the tree with the solvent code. No further
calibration of the retired prior is planned.

### Stage 6.5 — competition-path retirement experiment (planned)

The experiment that would conclusively retire the gridding + post-hoc NU
competition path (including `nu_refine` shell extension) in favor of
`rec_backend=pcg` with the `Q_NU` prior. Code preconditions are complete as
of 2026-08-29: every workflow (abinitio3D, refine3D, refine3D_auto) runs the
pcg arm with zero prior flags, no competition machinery executing anywhere
(startup, per-iteration, or postprocess), envelope masking available under
`automsk=yes`, and cross-backend FSC/cFAR comparability restored by the mask
unification. Experiment design and acceptance criteria: TO BE DEFINED by the
user before execution -- candidate axes are the Gate C fixture set plus
real-data refine3D_auto pairs (gridding+NU competition vs pcg+Q_NU) judged
at matched convergence on resolution, cFAR, map quality, over-regularization
diagnostics, and wall-clock/cost. On acceptance: remove the retired solvent
prior code, the post-hoc NU competition path for pcg, and revisit
`nu_refine`'s existence; record removals here.

### Stage 6.6 — always-on adaptive band granularity (IMPLEMENTED 2026-08-29; build/runtime gates user-side pending)

The replay-setting successor to `nu_refine`-style resolution extension,
designed and implemented 2026-08-29. The problem: the
solver's finest band spans 8 A to Nyquist as ONE weight block (the
evidence probes stop at ~5-4 A), so high-resolution refinement pays a
granularity tax -- one confidence governs the whole fine tail, and the
matching-lp handoff is quantized to the candidate ladder near the
frontier. The competition path solved this with the shell challenger
(`nu_refine=yes`), which had to be OPT-IN because filter-side acceptance
was irreversible within the iteration: a wrongly accepted shell passed
its noise directly into the references.

Design: make the band structure a deterministic function of the current
evidence state, ALWAYS ON -- no flag, no successor to `nu_refine=yes|no`.
Keep the static 20/12/8/5 A ladder and extend it geometrically (the
existing ~0.6-0.67 per-step spacing) until the next boundary would pass
max(current base FSC=0.143 crossing x margin, crop Nyquist), with
matching evidence candidates per new boundary. Band count then grows with
resolution the way the abinitio crops do; the structure is derived from
the pair, shared by both halves, and frozen before either replay, per the
existing evidence contract.

Why always-on is safe where the shell challenger was not (risk
inversion): subdivision has NO acceptance step. A new finest band earning
low confidence adds penalty to the fine tail (the conservative direction
-- today's one-block band 4 under-penalizes fine noise inside a
partially-evidenced region); earning high confidence, the penalty
retreats exactly where evidence supports it. Nothing is passed, nothing
is cut, and the structure is refrozen from fresh evidence every
iteration, so a poorly placed boundary self-corrects.

Implementation surface is evidence-side only: the solver
(`set_nu_prior`, band index, operator, preconditioner approximation) is
already nbands-generic. Changes: derive the ladder/limits per iteration
(replacing the `NU_EVIDENCE_BAND_LIMITS` constant as the sole source),
pack them in the compact state, generalize
`expand_nu_evidence_band_weights` past the fixed four, and thread the
state-carried limits to `set_nu_prior` (the provenance string already
records `bands_A`/`candidates_A`, so the evidence identity hash tracks
the structure for free). Side benefit: the matching-lp handoff steps
refine with the ladder (the 4.01/4.12/4.42 A staircase in the msp1 log
smooths out).

Gates: the Gate A partition test at each band-count transition (the bank
is not a tight frame; refining the partition changes the operator --
measured, not assumed); the suppression readout watched across
transitions (a mid-refinement subdivision redistributes where
suppression concentrates); one lambda-ladder recheck at high band count.
Cost note: two padded FFT pairs per extra band per application (~75%
more Q_NU application cost at ~7 bands), fixed per solve and independent
of particle count; a few extra candidates in the once-per-iteration
evidence analysis.

IMPLEMENTATION RECORD (2026-08-29). As designed, evidence-side only:

- Constants `NU_EVIDENCE_MAX_NBANDS=8`, `NU_EVIDENCE_BAND_RATIO=0.64`,
  `NU_EVIDENCE_FRONTIER_MARGIN=0.8` and the module-state active ladder
  `nu_evidence_band_limits_active` in `simple_nu_filter`;
  `nu_evidence_summary%band_limits`/`supported_fraction` are now
  allocatable (band count is part of the frozen state).
- `derive_adaptive_evidence_bands` (bank submodule): always sets the
  active ladder -- static without a frontier (exact pre-6.6 behavior,
  which is what every non-replay caller and `simple_test_nu_envmask`
  get), extended geometrically toward
  max(0.8 x frontier, 2 x smpd) otherwise, two candidate probes per
  appended band (geometric midpoint + boundary). `setup_nu_dmats` gains
  optional `evidence_frontier_lp` (requires `evidence_source`);
  `init_nu_filter` gains optional `extra_candidate_lps`, appended
  strictly monotone past the static ladder, grid-snapped duplicates
  dropped, Nyquist- and cap-clamped.
- `build_nu_evidence_state` sizes band support, summary, provenance
  (`bands_A` now dynamic), and the identity checksum from the active
  ladder; `nu_evidence_state_is_valid` checks band count in
  [4, NU_EVIDENCE_MAX_NBANDS], strictly decreasing positive limits, and
  that the first four entries ARE the static ladder;
  `expand_nu_evidence_band_weights` and the solver
  (`set_nu_prior`/band index/operator/preconditioner) were already
  nbands-generic and are untouched.
- `build_nu_replay_evidence` takes the evidence pair's FSC=0.143
  crossing as the frontier and returns the state-carried `band_limits`
  alongside `band_w`; both PCG execution paths thread
  `res0143s(state)` in and pass the returned ladder to `set_nu_prior`,
  so operator and evidence can never disagree
  (`NU_EVIDENCE_BAND_LIMITS` no longer reaches the strategy). The log
  line `>>> NU ADAPTIVE EVIDENCE BANDS: <n> bands, finest boundary ...`
  fires whenever the ladder extends.
- No new parameter anywhere; `nu_refine` has no successor knob.

Runtime observables for the user-side gates: the adaptive-bands log line
appearing once resolution passes ~6 A (finest static boundary 5 A,
margin 0.8: first extension at frontier ~4 A / floor permitting);
`pcg_nu_supported_fraction_band05+` lines and a longer `bands_A` list in
the evidence block; smoother matching-lp handoff steps; suppression
percentage continuity across a band-count transition.

### Stage 6.6 first run: 1WCM truth-judged (2026-08-29) — unguarded subdivision DEGRADES; evidence-gated retention added

First adaptive-bands run, same fixture and budget as the Stage 6 Gate C
runs. The mechanism fired exactly per the rule: frontier 2.943 A ->
5 bands, finest boundary 3.2 A (next step 2.05 A below the 2.35 A floor,
correctly stopped); 10 candidates (the appended midpoint probe snapped
to the existing 4 A find and was correctly deduplicated, only 3.2 A
appended); dynamic `bands_A`/`candidates_A` in the provenance; all
harness gates PASSED. But the truth verdict is a clean DEGRADATION
against the 4-band control across the entire fine band, dropping below
even the gridding reference from ~3.7 A (truth-FSC at k=79/89/99/109:
0.739/0.440/0.219/0.098 vs the 4-band 0.813/0.550/0.308/0.148 and
gridding 0.766/0.472/0.227/0.117).

Diagnosis, visible in one line: `band05 support = 0.000000`. In the
solver partition the new band 5 spans EVERYTHING finer than 5 A, so the
zero-support subdivision replaced the fine tail's partially evidenced
weight ((1-0.19)^2 ~ 0.65) with the full penalty 1.0 -- over-suppressing
genuine 3-5 A signal band-wide. The design argument that low confidence
in a new band is "the conservative direction" was right about direction
and wrong about magnitude. The zero reading also exposes a detector
limitation: band confidence is softmax WINNER MASS on candidates at
least that fine, and at fine candidate spacing (3.2 vs 4.0 A) the
whitened objective barely separates neighbors, so the mass dilutes
toward zero at the frontier even where signal demonstrably extends (the
pair's crossing is 2.94 A). Suppression readout moved 35.9 -> 31.3%
across the partition change -- the predicted discontinuity, now
measured.

Fix (implemented): EVIDENCE-GATED RETENTION, the always-on form of
"subdivide when evidence supports it". After band support is computed,
appended bands are pruned finest-first unless their mean support reaches
`NU_EVIDENCE_MIN_BAND_SUPPORT = 0.01`; the static four are never pruned,
so pre-6.6 behavior is the guaranteed floor and a zero-support
subdivision self-neutralizes to the 4-band operator (logged as
'pruned to N band(s)'). Still no flag, no acceptance step on KEPT bands,
per-iteration refreeze unchanged.

VERIFIED (2026-08-29, same-day rerun): the guard behaves exactly as
designed. Extension proposed (5 bands, 3.2 A), then
'pruned to 4 band(s); appended band(s) earned no support'; `bands_A`
back to the static four while `candidates_A` retains the 3.2 A probe
(which keeps informing selected cutoffs and the LP handoff without
penalizing); suppression restored to 36.0/36.9% (4-band control:
35.9/36.8 -- continuity recovered); shipped-pair crossings identical to
the 4-band control (3.241/2.844); truth-FSC matches the 4-band control
to the fourth decimal across the band; all gates PASS. The retained
extra candidate perturbs the softmax slightly (band01 support 0.9227 vs
0.9279, temperature 4.18e-3 vs 6.04e-3) with no measurable truth effect
-- the acknowledged and acceptable footprint of a pruned proposal.

OPEN (Stage 6.6b, if adaptive granularity is to ever engage usefully):
the fine-band confidence criterion needs a spacing-aware form -- winner
mass systematically underestimates support as candidate spacing shrinks;
a margin-based criterion (as the envelope uses) or cumulative-mass
calibration against the null are the candidates. Until then, adaptive
banding is expected to read as a no-op on most data (extension proposed,
then pruned), which is safe and honest.

### Stage 6.6 FINAL FORM (2026-08-29, user-directed): nu_refine mirror of the gridding challenger

The 6.6b spacing-aware criterion (a pairwise-margin sigmoid) was briefly
implemented and WITHDRAWN unmeasured on user direction: too complicated,
and unnecessary -- the codebase already contains the proven answer. The
final design mirrors the gridding path exactly:

- The a-priori frontier-tracked candidate proposal (the
  `evidence_frontier_lp` plumbing) is REMOVED. Candidate-bank extension
  in the evidence analysis is instead the existing nu_refine
  high-resolution shell walk (`extend_nu_filter_highres_shells`): one
  Fourier shell at a time from the populated frontier, accepted only on
  strict unary WIN-FRACTION at the tested frontier (>= 5% plus an
  absolute seed floor) -- the spacing-robust criterion that has gated
  gridding resolution extension all along, and precisely why winner-mass
  confidence is sound here: a band is only ever created over candidates
  that WIN somewhere, and dilution only zeroes dominated candidates.
- Gated by the EXISTING `nu_refine` flag, with its established meaning
  and workflow split: abinitio3D keeps the discrete static ladder (its
  stage policy pins nu_refine=no; the validated default-on Q_NU
  configuration is untouched); refine3D_auto restores its unconditional
  nu_refine=yes default on BOTH backends -- the gridding filter
  challenger there, the Q_NU evidence-bank extension here. No new
  parameter; the earlier backend-conditional nu_refine default is
  reverted.
- The walk runs inside `build_nu_replay_evidence` between
  `optimize_nu_cutoff_finds` and `build_nu_evidence_state` (the exact
  sequence of refine3D_auto's gridding bootstrap), logged as
  'EVIDENCE BANK EXTENDED BY n ACCEPTED SHELL STEP(S) TO x A'.
- The band ladder is derived from the ACTUAL bank at evidence-build
  time: the static four bands, extended geometrically (x0.64, cap 8)
  only while an accepted candidate is at least as fine as the next
  boundary. Static bank => exactly the four-band behavior; no band can
  exist without a challenger-validated probe. The evidence-gated
  retention guard stays as belt-and-braces behind the challenger gate.
- Guards updated: `nu_refine=yes` on the pcg backend now requires the
  Q_NU replay (clear hard error otherwise, replacing the blanket
  'does not yet support' errors); the post-hoc P_tau+NU pcg path keeps
  no extension support.

VALIDATED on the 1WCM harness (2026-08-29, both runs, all gates PASS):

- nu_refine=no control: EXACT parity with the verified four-band Q_NU
  Gate C result -- support fractions, suppression (35.9/36.8%), shipped
  crossings (3.241/2.844), and truth-FSC identical to all printed
  digits. With no walk there is no candidate footprint at all, so
  abinitio3D's discrete-ladder mode is the validated configuration
  verbatim.
- nu_refine=yes: the walk accepted 30 shell steps to 2.72 A (bank at the
  24-candidate cap), band 5 (3.2 A boundary) earned REAL support (0.064,
  vs 0.000 under the a-priori proposal) because it sits over
  challenger-validated candidates that win at frontier voxels; ladder
  correctly stopped at 5 bands (next boundary 2.05 A below the finest
  accepted candidate). Truth verdict: IMPROVEMENT over the already
  dominant Q_NU control at every shell finer than ~3.3 A (2.91 A: 0.587
  vs 0.578; 2.61 A: 0.344 vs 0.335; 2.49 A: 0.237 vs 0.227; gains
  persist to the band edge) against a ~0.003 truth-FSC dip in the
  3.8-4.4 A region -- the intended granularity trade, favorable for a
  refinement tool. Suppression continuous (34.6/35.6%); shipped-pair
  inflation 4.6% vs the control's 3.4% (diagnostic only, inside
  guidance).
- One defect found and fixed by the first nu_refine=yes attempt: the
  provenance and identity-seed buffers (LONGSTRLEN=1024) overflowed with
  the 24-candidate provenance ('End of record' at the identity internal
  write); both are now XLONGSTRLEN.

### Run record (2026-08-29): PfCRT abinitio3D -- Q_NU INERT at lambda_rel=0.1 on real data

First real-data abinitio3D with the Q_NU replay on PfCRT (EMPIAR-10330
final refined particle stack, mskdiam 160 A;
`Processing/pcg_integration/PfCRT/ABINITIO3D_OUTPUT_RESTART1`, code at
the pre-6.6 state -- equivalent here since abinitio3D pins
nu_refine=no). Outcome: map quality NOT competitive with gridding + NU
filter. Diagnosis from the output doc:

- Stages 1-5 (LP 20 -> 11): Q_NU off by design (filter mode none),
  mirroring the gridding path's NU-off early stages. Correct behavior.
- Stages 6-8 (LP 9.7 -> 6.0): Q_NU on at the dynamic default
  lambda_rel=0.1, evidence rebuilt per iteration with sensible band
  trajectories (band 1-3 support growing 0.76->0.83 / 0.69->0.78 /
  0.15->0.67; band 4 low, 0.001->0.05), BUT prior-energy suppression
  never exceeded 4.8% (0.3% at NU onset, 1.6% at the final rec) vs ~35%
  on the 1WCM harness and the bgal run at the same lambda_rel. The
  `PCG NU PRIOR INERT (< 5%)` diagnostic fired at essentially every
  NU-era iteration, recommending ~3x lambda.
- Failure mechanism: on the pcg path Q_NU is the ONLY spatial
  regularization (post-hoc NU filtering skipped by design), so an inert
  prior means refinement aligns against effectively unshaped ML maps
  while the gridding path aligns against NU-filtered references. The
  data-scale anchoring of lambda_eff evidently does not normalize
  regularization strength across datasets.
- Secondary observations: null calibration much noisier than 1WCM
  (null_bias_median 0.46 vs 0.06), only 9 candidates survived, and the
  matching low-pass sat pinned at the static bank's 3.98 A floor
  through stage 8.

Next step (user-directed): rerun PfCRT at higher lambda_rel (0.3, then
1.0 if suppression still lands under ~15-20%; suppression is sublinear
in lambda so the 3x hint may be conservative) to establish empirically
what a suitable suppression target is. A suppression-targeted
closed-loop auto-lambda (the replay already measures suppression every
reconstruction) is the candidate long-term fix, to be designed only
after the scan fixes the target band.

### Q_NU operator performance pass (2026-08-30, operator-preserving)

`apply_nu_precision` cost 3*nb+1 padded FFTs per matvec (16 at nb=5)
against the kernelized data operator's 2 -- the prior was ~8x the data
term per PCG iteration. Two operator-preserving fixes landed (identical
operator up to floating-point reordering):

- **Fourier-side adjoint accumulation:** the bands are disjoint radial
  partitions, so the per-band masked spectra are assigned into one
  accumulator and synthesized with a SINGLE final inverse transform
  (IFFT is linear). Saves nb-1 inverse FFTs.
- **Shared forward FFT (kernel path):** the kernel matvec stashes its
  pristine FFT(pad(p)) via `stash_nu_forward` before the Khat multiply
  (deapodization on, so the spectrum is of the iterate itself); the
  prior completes the mean-centering in Fourier space with the
  precomputed box-window spectrum `nu_winhat = FFT(pad(1))`, since
  FFT(pad(p - m)) = FFT(pad(p)) - m*FFT(pad(1)). Saves the prior's own
  forward FFT. The matrix-free path cannot share (its spectrum is of
  the deapodized iterate) and keeps its own forward.

Net: 3nb+1 -> 2nb+1 padded transforms on the kernel path (16 -> 11 at
nb=5), 2nb+2 on matrix-free. Allocator churn also removed: persistent
workspaces (`nu_cmat0/nu_cacc/nu_cband/nu_xb`, `ensure_nu_workspaces`)
replace the per-matvec allocations and all hot-path `pad_vol`/`crop_vol`
temporaries (direct `get_rmat_ptr` interior writes on a zeroed wimg).
Back-pocket option NOT taken (operator-changing, needs its own 1WCM
validation): running the radial band frame at the native box instead of
the padded lattice (~8x cheaper FFTs at padf=2). Memory flag stands:
`nu_band_w` is box^3 x nb (~2 GB at box 400, nb 8).

Verification requirement (R5-adjacent): rerun the 1WCM Q_NU control --
suppression, support fractions, and truth-FSC must reproduce the
verified values to printed digits (FP reordering may wiggle the last
digit); any larger deviation means the rewrite changed the operator.

**PARITY VERIFIED (2026-08-30, 1WCM control, all gates PASS):**
suppression 35.89/36.82 and shipped crossings 3.241/2.844 identical to
the verified control to all printed digits; all 129 truth-FSC shells
bit-identical as printed; support fractions differ only in the sixth
decimal (0.927938 vs 0.927939 etc.) -- exactly the FP-reordering
signature expected from Fourier-side accumulation. Wall-clock not
comparable on this run (host under load ~11 from concurrent jobs);
judge the speedup from NU-era iteration timings in the next quiet
PfCRT/refine3D run.

### Active dev list (2026-08-29, user-prioritized)

Follow-on note: some PfCRT maps look acceptable, but with the prior
inert that is the unregularized ML solve doing the work, not Q_NU --
the acceptable-looking outputs do not validate the prior.

1. **PRESSING -- NU-evidence locscale-style nonuniform postprocessing**
   (`nu_evidence_local_sharpening.md`). A single isotropic B-factor
   does not work for most maps other than bgal and streptavidin; graded
   local amplitude scaling from the frozen NU evidence is the
   model-free fix. The proposal's Gate C/D precondition is now met
   (Q_NU cleared its gate program and the Stage 6.6 walk validated), so
   implementation is unblocked and elevated to the top of the queue.
   **IMPLEMENTED (2026-08-29)** as the isolated `postprocess_nu`
   commander + `simple_nu_filter_sharpen` submodule (user-directed
   isolation from the standard postprocess path). v1 (confidence-gain
   band stack) over-sharpened the PfCRT core and is RETIRED as a
   recorded failure; v2 is the classical shrink-then-sharpen localized
   by the evidence (Guinier B inside the evidenced local Butterworth
   passband, null-flattened solvent). Design and records in
   `nu_evidence_local_sharpening.md` §3b/§3c; validation per its §4
   plan pending.
2. **Auto-lambda design after regularization parameter testing.** The
   PfCRT lambda_rel scan (0.3, then 1.0; in flight) establishes the
   target suppression band empirically; only then design the
   suppression-targeted closed-loop lambda_eff controller (the replay
   already measures suppression every reconstruction). R9: record the
   target band here before implementing the controller.

   Scan record and root-cause analysis (2026-08-30). PfCRT suppression:
   lambda_rel=0.1 -> 3.25%, lambda_rel=0.3 -> 8.6% (shipped-pair 3.93 A,
   diagnostic banner NOMINAL). The x2.65 response to x3 lambda is almost
   perfectly LINEAR -- the weak-prior perturbative regime, so 1WCM-like
   shaping (~35%) needs an order of magnitude, not small multiples;
   next points lambda_rel=1.0 and 3.0 to bracket, and where the
   response bends sublinear is the saturation information the
   controller needs. Root cause of the cross-dataset failure:
   `update_lambda_from_density` anchors `data_scale` to the LOW-BAND
   mean of the raw data diagonal D, so lambda_eff = lambda_rel *
   mean(D_lowband) -- but the prior competes against D in the FINE
   shells where the unsupported-band content lives. The effective
   strength is lambda_rel * D_low/D_fine, and that spectral falloff
   ratio is a dataset property (CTF envelopes, particle count, ice,
   crop grid): gentle on the 1WCM phantom (0.1 -> 35%), ~10x steeper
   on PfCRT (0.1 -> 3%). Candidate fixes, to be decided AFTER the
   scan: (a) anchor data_scale to the prior-active shells (mean D over
   shells finer than the coarsest band boundary, or D weighted by the
   prior's spectral footprint), making lambda_rel dataset-invariant by
   construction; (b) close the loop on measured suppression.
   Discriminating test once a ~30-40% PfCRT point exists: does
   fine-shell anchoring alone collapse the PfCRT and 1WCM operating
   points onto one lambda_rel?

   **TARGET RECORDED (2026-08-30, user, R9): aim for ~60% prior energy
   suppressed.** Empirical basis: at that operating point the NU
   effects are clearly visible on real data (msp1). This sits at the
   top of the diagnostic banner's NOMINAL band (5-60%), so the banner
   thresholds should be re-centered around the 60% target when the
   auto-lambda controller is designed; the controller's setpoint is
   60%, not the 1WCM harness's ~35%.

   **AUTO-LAMBDA IMPLEMENTED (2026-08-30).** Design recorded before
   first run (R9):
   - Plant model: one-pole amplitude-suppression law
     s = g*lambda/(1+g*lambda), validated on the scan (PfCRT g=0.336
     from the 0.1 point predicts s(0.3)=9.2% vs 8.6% measured; 1WCM
     g=5.4 -- a 16x cross-dataset gain spread, which is why no fixed
     lambda_rel can work).
   - Controller: memoryless one-step secant per refinement iteration,
     `resolve_nu_autolambda` in `simple_rec3D_pcg_strategy` (both exec
     paths, before attachment): reads the previous iteration's
     `PCG_NU_STATS_FILE` (lambda used + state-mean suppression, the
     file the replay already persists), identifies g, solves for the
     60% target. Deadband hold at 60 +/- 5%; multiplicative step clamp
     x5; lambda_rel bounds [0.01, 30]; suppression floored at 0.1% and
     capped at 99% for finite identification. Both halves and all
     states share one lambda, same freezing discipline as the evidence.
   - Activation: ONLY when pcg_nu_lambda_rel was left to its dynamic
     default (`l_pcg_nu_autolambda`, set in the parameters dynamic
     default block). An explicit strength pins lambda -- every recorded
     control, harness run, and ladder point stays reproducible, and
     explicit 0 = P_tau control is untouched (R5/R10 preserved). Cold
     start (no stats file) keeps the 0.1 default, so single-shot
     reconstructions in fresh directories are deterministic.
   - Diagnostics: stats file gains a PCG_NU_AUTOLAMBDA key; the
     convergence banner is re-centered on the target (INERT <5%, BELOW
     TARGET <45%, ON TARGET 45-75% aim 60%, OVER >75%) and, in auto
     mode, reports the lambda in use instead of advising manual
     changes; the controller logs MEASURED/HOLDING/adapted-lambda lines
     per reconstruction.
   - Deliberately NOT done: fine-shell re-anchoring of data_scale
     (operator-meaning change, would force re-baselining per R2; the
     closed loop absorbs the anchoring error in one update). Recorded
     as the cold-start improvement to revisit only if convergence to
     target proves too slow on some dataset class.
   - Validation plan: (i) 1WCM harness with explicit lambda -- must be
     bit-identical to the verified control (controller inert when
     pinned); (ii) PfCRT/msp1 refine3D or abinitio3D WITHOUT a lambda
     flag -- expect convergence to 55-65% suppression within ~3
     iterations from cold start, then station-keeping; watch the
     shipped-pair inflation diagnostic while the controller holds.

   **FIXED 60% SETPOINT FALSIFIED -- PfCRT REGRESSION (2026-08-31).**
   First production abinitio3D runs of the auto-lambda controller
   (c2b419ad, no lambda flag, nu_refine=yes): PfCRT regressed from
   3.93/3.98 A shipped pairs (b082ef4a, pinned lambda_rel=0.3, supp
   8.6%) to 5.86 A with the transmembrane helices resolved as sausages
   -- the controller drove lambda toward the 60% setpoint (g=0.336
   implies lambda ~4.5, ~15x the proven strength). streptavidin, msp1
   and FlhB held or improved at the same commits. Conclusion: the
   plant GAIN spread (16x) that motivated auto-lambda extends to the
   TARGET -- the good operating point is a dataset property (PfCRT
   ~9%, 1WCM ~35%, msp1 ~60%), so no fixed suppression setpoint
   transfers.

   **AUTO-TARGET IMPLEMENTED (2026-08-31, user-directed).** The
   setpoint itself is now controlled; design recorded before first run
   (R9):
   - Structure: cascade. Inner loop unchanged (one-step secant
     auto-lambda tracking the setpoint). New outer loop
     (`resolve_nu_supp_target` logic inside `resolve_nu_autolambda`)
     owns the setpoint, driven by the shipped-pair FSC=0.143
     trajectory -- the persisted over-regularization diagnostic, never
     a resolution claim.
   - Law: AIMD (additive-increase / multiplicative-decrease). Per
     reconstruction, comparing the previous iteration's state-mean
     shipped-pair crossing to the one before it (2% relative
     deadband): improved -> setpoint +5 percentage points; degraded ->
     setpoint x0.6; stalled -> hold. lp-limited stage plateaus stall
     and therefore hold; crop/box changes that perturb the crossing at
     stage transitions can only trigger a conservative back-off.
   - Bounds [5, 75]% (banner inert floor to banner over ceiling); cold
     start 15% -- gentle by construction (above PfCRT's proven 8.6%,
     ramping toward msp1-class operating points only while the shipped
     pair keeps improving: 60% is reachable in 9 improving
     iterations).
   - State: memoryless via `PCG_NU_STATS_FILE`, which now also carries
     PCG_NU_SUPP_TARGET, PCG_NU_AUTOTARGET, PCG_NU_SHIP0143_AVG and
     the one-step history PCG_NU_SHIP0143_PREV (lifted before each
     rewrite).
   - Activation: ONLY when `pcg_nu_supp_target` is left to its dynamic
     default (`l_pcg_nu_autotarget`); an explicit value in [5,75] pins
     the setpoint (reproducible controls, R2), and requires the
     auto-lambda controller (an explicit lambda pins the strength
     outright, hard error otherwise). Registered on refine3D and
     abinitio3D (forwarded to the refine3D stages; stripped with the
     other PCG-only keys off non-PCG child command lines).
   - Diagnostics: the convergence banner reports AUTO-TARGET/PINNED
     TARGET with the live setpoint; the manual-lambda advisory bands
     are re-centered on the reported setpoint (default 60 when
     absent); the controller logs IMPROVED/DEGRADED/STALLED setpoint
     lines per reconstruction.
   - Validation plan: (i) pinned lambda -- both controllers inert,
     bit-identical to the verified control; (ii) PfCRT abinitio3D with
     no flags -- expect the setpoint to hold near the cold start while
     shipped pairs march, back off on any inflation, and the run to
     recover the ~3.9 A pinned-lambda result; (iii) msp1 -- expect the
     setpoint to ratchet toward its known-good ~60% while shipped
     pairs improve, with no regression vs the fixed-60% runs.

   **FINAL-RECONSTRUCTION Q_NU POLICY (2026-08-31, user-directed).**
   The original-sampling final reconstructions of abinitio3D and
   refine3D_auto previously dropped the PCG backend and its prior
   entirely (abinitio3D: `prep_final_rec_cline` rebuilt the child line
   without `rec_backend`, shipping a gridding P_tau final map;
   refine3D_auto: a plain unregularized objfun=cc rec3D). The
   objfun=cc passes exist for sigma availability -- registration-box
   sigmas are crop-incompatible at the final box, which is exactly
   what `bootstrap_rec3D` solves (cc pass -> half-map sigma
   derivation -> euclid ML pass). Policy now: **the final map is
   reconstructed on the refinement's backend, and on pcg the Q_NU
   replay regularizes it in-solve.** Implementation:
   - `bootstrap_rec3D` strips the Q_NU keys off its unregularized
     sigma pass (no ml_reg, activation contract) and lets them and
     filt_mode flow into the regularized pass.
   - abinitio3D `prep_final_rec_cline` forwards
     rec_backend/maxits_pcg/rtol and, when the final stage uses
     ml_reg, sets filt_mode=nonuniform and forwards the pinning keys.
   - refine3D_auto's final rec now runs `bootstrap_rec3D` (previously
     unregularized cc), with bootstrap sigmas written at endit+2 so no
     crop-box sigma star is overwritten, and the same pcg/Q_NU
     forwarding.
   - Grid-transition calibration (corrected 2026-09-02): `lambda_rel`
     does not transfer unchanged across the staged-crop -> native-grid
     transition. The change in downscaling changes the effective Q_NU
     plant gain even though the suppression target remains meaningful.
     With auto-lambda active, `bootstrap_rec3D` therefore retains the
     learned suppression target, lets the unregularized sigma pass clear
     the old-grid response, runs one current-grid Q_NU calibration solve
     without postprocessing, and then reruns the regularized reconstruction
     after `resolve_nu_autolambda` has adapted from that measurement. Only
     the calibrated replay is shipped. A missing target starts at 15%.
     Explicit lambda controls skip calibration and remain pinned as usual.

   **PfCRT REGRESSION ROOT CAUSE -- HANDOFF GATE + FROZEN-EVIDENCE
   DEADLOCK (2026-08-31).** The post-auto-target PfCRT abinitio3D run
   (9b2424b5) still stalled at 5.86 A with the matching low-pass pinned
   at 5.97 A, "RIDING FROZEN EVIDENCE (no resolution advance)" and
   "SHIPPED PAIR STALLED" every iteration -- the controllers were
   holding on a stall they could not break, so the prior strength was
   never the driver. Git bisection of the record: both offending
   changes landed together in a6a5cc1e/9f3c3c9e ("PCG NU prior
   optimizations for speed"), the first commit after the good pinned
   b082ef4a runs (3.93/3.98 A):
   - The LP-set matching handoff switched from the raw finest selected
     per-voxel cutoff (`minval(cutoffs)`) to the 5% assignment-support
     percentile (`nu_evidence_finest_supported_lp` at
     `NU_ALIGN_LP_MIN_ASSIGNED_PCT`). On a small membrane protein only
     a small core carries fine evidence while the micelle belt
     dominates the assigned support, so the gate collapsed the matching
     bandwidth to ~the FSC crossing. streptavidin/msp1/FlhB resolve
     compactly and uniformly, so the percentile ~= the minimum and they
     were unaffected.
   - The frozen evidence cache rebuilds on FSC=0.143 ADVANCE -- but the
     alignment search is capped at the evidence-derived low-pass, so no
     advance can occur: a self-sealing loop (the age-5 forced rebuild
     regenerates from maps aligned under the cap and the gate
     re-collapses the handoff). The plateau then reads as convergence.
   Fixes (both, user-directed):
   - Handoff decoupled from the support gate: raw finest selected
     cutoff (`min_pct=0`), restoring the proven b082ef4a behavior.
     Sparse-but-real fine evidence must be allowed to pull the search
     band forward; over-extension widens the search, over-restriction
     deadlocks it. The gridding-path gate
     (`get_nu_filtmap_finest_selected_lp`) is untouched.
   - Binding-band rebuild condition in `nu_evidence_needs_rebuild`: the
     cache entry records the Fourier index of the matching low-pass it
     handed off (`handoff_find`); once the FSC crossing comes within
     one shell of that cap, the alignment band is the binding
     constraint and the evidence rebuilds from the live pair (with the
     nu_refine shell walk re-attempted). Principle: the cache may trade
     staleness for speed in the Q_NU band weights, but the SEARCH
     BANDWIDTH must never be governed by a frozen statistic.
   - Validation plan: PfCRT abinitio3D, no prior flags -- expect the
     matching low-pass to extend ahead of the FSC again and the run to
     recover the ~3.9 A result; streptavidin/msp1/FlhB reruns -- expect
     no change (gate inert for compact particles, cache rides only
     while non-binding, so the speed win survives where it was ever
     legitimate).

   **MATCHING-REFERENCE REGRESSION -- PLAIN NONUNIFORM TOPOLOGY
   (2026-08-31).** First refine3D_auto+pcg PfCRT run (auto-lambda and
   auto-target behaving: supp 13.3%, lambda 0.177, target held at the
   15% cold start): shipped pair stuck at the gridding run's
   ITERATION-1 state (FSC=0.500 at 7.96 A / 0.143 at 4.14 vs gridding
   3.98/3.61 by iteration 2), postprocessed maps massively
   over-sharpened. Root cause chain:
   - The matcher's plain-nonuniform contract consumes INDEPENDENT
     even/odd `_nu_filt` references with NO further filtering
     ("filtering done when volumes are assembled",
     `simple_matcher_refvol_utils`). The redundancy policy
     (2026-08-28: Q_NU in-solve => no post-hoc NU filter, no _nu_filt
     products) silently turned the matcher's first-iteration raw-refs
     fallback into the PERMANENT state: every iteration matched
     independent raw Q_NU half maps, unfiltered, at the evidenced
     matching lp (~4 A) from iteration 1.
   - abinitio3D never sees this because GOLD_STD_STAGE=TURNED_OFF
     rewrites every nonuniform stage to nonuniform_lpset -- MERGED
     single-reference matching, no independent per-half registration.
     The redundancy policy was only ever validated in that topology.
   - Independent per-half matching against each half's own
     unsuppressed in-band noise (mild auto-targeted replay: 13%
     suppression barely regularizes the refs; the old fixed-60%
     setpoint was inadvertently cleaning them) overfits each half to
     its own noise -> mid-resolution half-map divergence with a
     matched-noise 0.143 tail. Textbook gold-standard overfitting.
   - Over-sharpening: abinitio3D pins bfac=0; refine3D_auto lets
     postprocess auto-estimate B by Guinier on the SHIPPED map. Q_NU
     (and ML) amplitude suppression steepens the Guinier slope ->
     auto-B strongly negative -> massive sharpening.
   Fixes (2026-08-31, both implemented):
   - `filter_pcg_nonuniform_maps`: the Q_NU skip-everything branch is
     now LP-SET-ONLY (merged-reference topology keeps handoff-only
     behavior). Plain nonuniform falls through and generates the
     derived `_nu_filt` matching references (base = _unfil pair under
     ml_reg; the ML/Q_NU pair enters as the finest-bank aux
     replacement only when nu_refine=no, per the aux-channel
     contract); the matching-lp handoff is then the filter-bank
     finest selected lp, as in gridding. Shipped primary outputs
     remain the raw Q_NU maps.
   - `postprocess_volume_from_files`: automatic B-factor is now
     always estimated from the `_even_unfil`/`_odd_unfil` pair
     average when present (prior-free amplitudes; the estimate
     targets the underlying signal decay), falling back to the
     shipped map only when no unfil pair exists.
3. **Test nu_refine=yes with rec_backend=pcg in abinitio3D.** The
   nu_refine=no rationale was gridding-specific (the ML-regularized aux
   competitor supplied beyond-bank resolution implicitly,
   `l_use_aux = l_ml_reg .and. .not. l_nu_refine`); the pcg path has no
   aux channel, so nu_refine=no is a hard 4 A evidence/matching-lp
   ceiling. The walk is gold-standard-gated so early noisy stages
   should accept nothing; verify on the bgal replay first, watching
   per-stage extension lines against the stage lp schedule.
4. **Solvent constraint via the lowest-bin convention (IMPLEMENTED
   2026-09-01, user-directed; design recorded before first run).**
   Context: the
   2026-08-05 fundamentals changes removed the envelope from the NU
   filter (spherical-only `setup_nu_dmats`, a31c7e8ad) and flipped
   refine3D_auto to spherical FSC (envfsc=no, c7dfd5055); the retired
   solvent prior left the PCG solve with no solvent constraint at all.
   envfsc=yes is restored as the refine3D_auto hard default
   (2026-09-01): the density-envelope automask + phase-randomization
   corrected FSC (`evaluate_halfmap_pair`) again steers matching lp,
   NU gating, and convergence. This item drafts the remaining half:
   reinstating the solvent constraint in the NU filter and the Q_NU
   prior WITHOUT reopening the reason the envelope was removed.

   Design principles:
   - ONE mask, one convention, both consumers. The conservative
     density automask (`image_msk%automask3D`: lp-smooth -> otsu
     binarize -> connected components -> grow -> cos_edge) at the
     conservative `envmsklp`, computed ON THE FLY from the current
     base average at assembly time, never written to file for this
     purpose. NOT the aggressive NU evidence envelope (amsklp-scale),
     which keeps its narrower jobs (detergent-omitting
     matching-reference masking, diagnostics) and is never used for
     the FSC or the solvent constraint.
   - POST-ASSIGNMENT CLAMP, not objective support. The spherical
     support of `setup_nu_dmats` is an invariant (envelope-restricted
     objectives broke the null statistics -- the reason for
     a31c7e8ad). The objective, whitening, and null estimation stay
     spherical and untouched. The envelope enters only AFTER
     optimization:
     (a) NU filter: voxels outside the envelope have their selected
         label overridden to the COARSEST candidate (lowest-resolution
         bin) before ordered-label Potts smoothing, which then owns
         the boundary. New optional envelope argument on the
         label-selection/`optimize_nu_cutoff_finds` handoff; absent
         argument = current behavior, bit-identical.
     (b) Q_NU prior: `expand_nu_evidence_band_weights` assigns voxels
         outside the same envelope to the coarsest band at the
         maximum lack-of-evidence weight, so the replay applies its
         strongest fine-shell suppression in solvent -- the in-solve
         solvent-flattening constraint the retired Q_s solvent prior
         was reaching for, expressed through the existing band
         machinery (no new prior term, R10 mode-exclusivity
         untouched).
   - Implementation (as built): module-level clamp state
     (`nu_solvent_lmask`/`nu_l_solvent_clamp`) with public
     `set_nu_solvent_envelope`/`clear_nu_solvent_envelope`; armed by the
     caller after `setup_nu_dmats`, cleared by `cleanup_nu_filter`, so
     absent = bit-identical current behavior for every other caller
     (nu_filt3D, bootstrap refs, tests). Filter clamp in
     `optimize_nu_cutoff_finds` (label -> coarsest candidate before
     Potts, logged as ">>> NU SOLVENT CLAMP: N support voxels...").
     Evidence clamp in `build_nu_evidence_state` (solvent -> null
     assignment, zero band support, zero uncertainty; skipped voxels
     bypass the softmax entirely, so the summary statistics --
     null_fraction, supported_fraction, uncertain_fraction -- reflect
     the constraint; provenance gains `solvent_clamp=density_envelope`,
     which flows into the identity hash, so a frozen state built with
     the clamp can never be mistaken for one without). Armed at three
     sites: `filter_pcg_nonuniform_maps` (plain-nonuniform pcg),
     gridding volassemble `setup_nonuniform_filter`, and
     `build_nu_replay_evidence` (Q_NU evidence, frozen with the state)
     -- each computing the density automask on the fly from the base or
     evidence pair average at `envmsklp`. Mask cost: one automask3D per
     state per assembly (same order as the envfsc path already pays).
   - Freezing discipline: the Q_NU clamp mask is part of the frozen
     evidence state (rebuilt with it), so operator and evidence cannot
     disagree across the cache lifetime; the search-bandwidth binding
     rule from the frozen-evidence deadlock fix is unaffected (the
     clamp only ever COARSENS solvent voxels, never the molecular
     region the handoff reads).
   - Validation plan: (i) A/B against the envfsc-restored build alone
     (A = envfsc=yes only, B = A + solvent clamp) on PfCRT
     refine3D_auto, both backends -- the clamp's contribution must be
     measured, not confounded with the FSC restoration; (ii) 1WCM
     harness with the clamp active: truth-FSC must not degrade
     (constraint is solvent-only by construction); (iii) suppression
     readout shift: expect the measured %-suppression to RISE at
     unchanged lambda (solvent now contributes), so watch the
     auto-lambda/auto-target interplay -- the setpoint controller must
     not compensate the clamp away by lowering lambda; if it does,
     restrict the suppression readout to the in-envelope region.

   **FIRST-RUN FINDINGS AND FIXES (2026-09-01, streptavidin
   refine3D_auto+pcg).** Shipped pair 3.05-3.12 A, controllers and
   clamp active; two defects surfaced and were fixed:
   - Matching low-pass pinned at the static bank finest (4.036 A)
     while the evidence extended to 3.1 A: the plain-nonuniform
     matching-reference pass built only the static ladder. Fixed:
     `filter_pcg_nonuniform_maps` now runs the proven shell walk
     (`extend_nu_filter_highres_shells`) after optimization when
     nu_refine=yes, so the matching handoff advances with the
     evidence (">>> NU MATCHING BANK EXTENDED..." line).
   - The predicted controller interaction FIRED: total-energy
     suppression read 40% at lambda 0.1 (solvent term), auto-lambda
     drove lambda to the 0.01 floor, molecular regularization
     collapsed to ~6%. Fixed as designed: the suppression readout is
     now restricted to the evidenced region -- `reconstructor_pcg`
     accumulates exact per-band energies y_b^T W_b y_b total and over
     `l_nu_evidenced` (any input band weight < 1) during
     `apply_nu_precision`; `get_nu_prior_stats` returns the evidenced
     penalty and `report_nu_solve_stats` feeds the controllers the
     evidenced-region suppression.
   - Output tidied: routine NU diagnostics (Potts sweeps, whitening
     profile, envelope ICM sweeps, envmask stats block, evidence
     calibration/provenance details, solvent-automask banner lines,
     ML warm-start notices, prior-energy detail) demoted behind
     NU_DEV_OUTPUT. Kept by default: the low-pass assignment table
     (the per-iteration NU fingerprint, restored on review), the
     solvent clamp count, a one-line evidence summary (id + null +
     band supports), one line per half with lambda_eff + evidenced
     suppression, and the controller lines.
   - Distributed-exec ergonomics (2026-09-01, user-directed, applies
     to ALL distributed SIMPLE applications): partition scripts no
     longer tee their output into the master terminal -- parts append
     to SIMPLE_SUBPROC_OUTPUT only, and the scheduler prints one
     summary line per completed phase (">>> <PRG>: N PART(S)
     COMPLETED IN X s", simple_qsys_ctrl). Part-level errors are
     found in SIMPLE_SUBPROC_OUTPUT, as the existing hard-error
     messages already direct.
   - Master-phase utilization (2026-09-01, user-directed): during the
     master-side PCG solve/NU/evidence phases the partition workers
     are idle; on LOCAL execution the distributed master raises its
     OpenMP thread count to min(ncores, nparts*nthr,
     PCG_MASTER_NTHR_CAP=32) for the reconstruction phase and
     restores nthr before matching (never on a cluster; capped
     because OpenMP scaling saturates at these box sizes). Ordinary
     master work uses that team; concurrent half solves split it and
     bake the half-team size into their private FFTW plans.
   - **Concurrent distributed even/odd solves -- implementation pending
     runtime validation** (2026-09-02): the reverted whole-lifecycle
     OpenMP sections are replaced by explicit prepare/solve/finalize
     half jobs. Preparation remains serial and owns support-mask
     memoization, FFTW planning, raw-accumulator I/O, trailing-chain
     writes, prior attachment, and warm-start construction. Only two
     fully prepared, half-owned `solve_accum` calls enter concurrent
     sections; validation, NU firing statistics, diagnostics/logging,
     image construction, and teardown are serial again. Each operator
     now owns an explicit half-team FFTW plan size, and nested OpenMP
     levels are enabled only around the solve pair then restored. PCG
     profiling uses stateless local timestamps rather than the shared
     legacy timer. This removes the deterministic code path behind the
     prior box-300 failure: spherical `set_mask -> mask3D_soft` no longer
     runs in an OpenMP region. The shared/direct route remains serial because its
     builder image batch is not half-owned. Remaining gate: measure
     peak resident memory and wall time at boxpd=600, and retain the
     implementation only if the two prepared workspaces fit with useful
     speedup in the final abinitio3D reconstruction.
   - NU work deduplication (2026-09-01, user-directed): the Q_NU
     evidence phase and the matching-reference generation ran the
     full setup + solvent clamp + optimization + shell walk twice per
     iteration on the same base pair. With nu_refine=yes (identical
     setups: same _unfil pair, no aux channel), single state, and
     source=base_unfil, the evidence phase now RETAINS its optimized
     setup (retain_nu_filter_setup / nu_filter_setup_is_retained in
     simple_nu_filter; cleanup_nu_filter always clears retention) and
     filter_pcg_nonuniform_maps consumes it directly -- skipping
     setup, clamp, extension and envmask regeneration -- then tears
     it down. Retention is keyed on the finest-lp output being
     requested (the signal that a matching volassemble follows), so
     standalone reconstructions never leak a retained setup.
     Multi-state and nu_refine=no keep the two-pass behavior (the
     module holds one setup; the aux channel makes the setups differ).
   - Second-run confirmation (streptavidin): matching bank extends
     (10/9 accepted steps to 3.1-3.2 A), matching lp advances via the
     support-gated promotion (3.43 A while the raw finest sits at 3.1
     A with ~1% assignment -- the gate working as designed), evidenced
     suppression keeps lambda off the floor spiral (settles ~0.01
     with supp ~10-13% against targets 9-14%), B-factor sane (-62 to
     -68). One controller wart fixed: the shipped-pair crossing is
     Fourier-shell-quantized and the 3.05<->3.12 A adjacent-shell
     flip re-triggered DEGRADED/IMPROVED steps; the AIMD comparison
     now stalls on any change of fewer than two shells.

5. **Direct PCG support constraint (EXPERIMENTAL, 2026-09-01,
   user-directed).** Distinct from every prior and from the solvent
   clamp: a real-space support constraint on the solve itself. The
   machinery pre-existed -- `set_mask` already solves the projected
   system (P H P) u = P b with the soft spherical support -- so the
   mode is: install an ARBITRARY [0,1] envelope as P.
   - `reconstructor_pcg%set_mask_volume(mskvol)`: caller-supplied
     volume as the support (validated, clipped to [0,1]); both solve
     entries now also project the INITIAL GUESS onto the support (a
     warm start's out-of-support content was previously carried into
     the output untouched -- latent even for the spherical case).
   - CLI: `pcg_mskfile=<vol.mrc>` on reconstruct3D and refine3D
     (pcg-gated, hard-errors without the pcg backend or a missing
     file). All PCG solve sites route through
     `set_pcg_solve_support`: pcg_mskfile when given, spherical
     mskdiam support otherwise -- workflow policies are untouched, so
     this is a pure play/experiment surface.
   - Composition: NOT a replay prior (R10 untouched); composes with
     P_tau/Q_NU (set_nu_prior already folds self%mask into the band
     weights, so the support automatically shapes the prior too), and
     with the trailing chain (constraint is solve-side only; the
     accumulators stay constraint-free, so masks may change across
     iterations and the mode switches off cleanly).
   - Caveats recorded: the naive half-map FSC is optimistic inside
     the support by construction -- judge with the phase-randomized
     envfsc; the mask is shared between halves (the classic solvent-
     flattening compromise); the Fourier-shell preconditioner does
     not commute with a real-space P (still valid, mildly
     suboptimal).
   - Deferred (design in the 2026-09-01 discussion): the soft-penalty
     variant (+ lambda_m diag((1-M)^2)) and the fixed-background
     focused mode (x = x_fix + P*delta, RHS = P(b - A x_fix)) -- both
     drop out of the same set_mask_volume machinery when wanted.
     Mode semantics for the record: HARD is a change of variables
     (x = P u), so out-of-support unknowns are REMOVED -- no strength
     parameter exists, conditioning improves by the support fraction,
     and a wrong mask deletes real density. SOFT is a Gaussian prior
     with precision lambda_m(1-M)^2, so out-of-support density is
     SHRUNK, not deleted, and survives where the data insists --
     graceful under mask error but reintroduces a dataset-dependent
     strength (the auto-lambda failure mode). The Q_NU solvent clamp
     is the band-limited middle: soft, fine-frequencies-only.

   **SOLVE SUPPORT POLICY (2026-09-01, user-directed; supersedes the
   first per-solve draft).** Principle: with PCG a mask belongs in
   the ESTIMATOR, not in post-processing -- post-hoc masking cannot
   improve a map, a support constraint can. Three rules (user,
   2026-09-01), over the two PCG passes (1) unfil and (2)
   regularized:
   - `automsk=yes` enables the conservative density envelope (Cyril's
     automask3D at envmsklp) as solve support at all.
   - `envfsc=yes` + automsk=yes: BOTH passes take the envelope.
   - `envfsc=no`  + automsk=yes: pass (1) keeps the SPHERICAL support
     (the FSC pair stays unconstrained), pass (2) takes the envelope.
   - Consequence, and the reason envfsc is the right switch: when
     pass (1) is envelope-constrained, the FSC pair is ALREADY
     masked by the estimator, so the envfsc masking +
     phase-randomization preprocessing in `evaluate_halfmap_pair`
     must not run again -- it is now gated on rec_backend (skipped
     for pcg, kept for gridding, where it remains the only way to get
     the envelope into the estimate). The envelope is still derived
     and returned there, because the automask artifact has other
     consumers (postprocess envfsc, matcher fallback, final rec).
   - ONE support envelope per state feeds both passes, so when both
     take it the FSC pair and the shipped pair stay on the same
     footing.
   - Derived from the reference volume the iteration matched against
     (lag-one, the same lag the matching references carry), which
     also resolves the ordering problem: the support exists before
     the base solves. Missing/startvol reference => spherical
     fallback.
   - The envelope is the CONSERVATIVE density automask at envmsklp
     (Cyril's 20 A mask), used as-is. An earlier draft coarsened it
     (max(30 A, 2*envmsklp)) out of concern for the NU evidence's
     null calibration and whitening; that was over-engineering
     (user, 2026-09-01): the 20 A envelope is already generous --
     protein plus micelle, dilated, soft-edged, and about half the
     spherical support on both datasets measured (PfCRT: 941k of
     1.93M support voxels outside it; streptavidin: 102k of 217k) --
     and NU refinement predates the shell whitening and works on
     plain Euclidean unaries, so the whitening MAD is not a design
     constraint. The evidence readiness contract (null_fraction
     outside [0.01, 0.90] hard-errors) is the guard if a specimen's
     envelope ever proves too tight for the null calibration; watch
     pcg_nu_null_fraction in the evidence banner on first runs.
   - Precedence in set_pcg_solve_support: explicit pcg_mskfile > the
     per-state density envelope > spherical mskdiam.
   - Backend split, for the record: gridding can only put the
     envelope into the MEASUREMENT (post-hoc mask + phase
     randomization); PCG puts it into the ESTIMATION (support
     constraint). Same envelope, same intent, different mechanism --
     and never both at once.

7. **refine3D_auto startup: reconstruct -> mask -> re-reconstruct ->
   refine (IMPLEMENTED 2026-09-01, user-directed).** Root cause of the
   PfCRT collapse, established by elimination (the envelope-masking and
   solver-convergence hypotheses were both falsified by controls --
   automsk=no reproduced the collapse exactly, and reconstruct3D at
   maxits_pcg=2 from the same orientations gives 5.97/3.93 A with cFAR
   0.78, better than refine3D_auto's first iteration): **the first
   matching had no proper references to match against.**
   - The matcher's nonuniform contract is that filtering is
     assembly-owned (`if( l_ml_reg .or. l_nonuniform_mode ) then !
     filtering done when volumes are assembled`). At iteration 1 the
     `_nu_filt` references do not exist yet, so it fell back to the RAW
     input half maps -- and then applied NO filter to them, matching
     unfiltered noise out to the FSC=0.143 band (~3.9 A on PfCRT, where
     that curve has already fallen to 0.16 and is about to cliff).
     abinitio3D's final stage, by contrast, matched a MERGED,
     NU-filtered reference at an explicit ~6 A lp.
   - The gridding path was immune because
     `prepare_nu_bootstrap_refs_from_raw_halves` generated NU-filtered
     startup references; it returned immediately for pcg ("the Q_NU
     replay regularizes in-solve"), which is true from iteration 2
     onward and false for iteration 1.
   - High-contrast specimens (bgal, streptavidin) survive matching
     unfiltered half maps; a low-contrast detergent-solubilized
     membrane protein does not. Hence "works on big things or things
     with no detergent".
   Implementation, three parts:
   (a) `reconstruct3D` on the pcg backend now produces the SAME
       reference products a refinement iteration does: `distr_execute`
       calls `filter_pcg_nonuniform_maps` after the master when
       filt_mode is nonuniform, generating the `_nu_filt` matching
       references and writing the evidence-derived matching-lp handoff
       into the project. Previously it called the master and returned,
       so no workflow could reconstruct-then-refine correctly.
   (b) refine3D_auto runs a STARTUP sequence before any matching:
       `calc_pspec` (particle-spectrum sigmas) then ONE regularized
       `reconstruct3D` carrying the refinement's own filt_mode/automsk/
       envfsc/nu_refine. It leaves the automask, the NU evidence
       envelope, the `_nu_filt` references and the lp handoff in place,
       and its output becomes vol1.
       `prepare_nu_bootstrap_refs_from_raw_halves` is deleted as
       superseded.
       SIGMA BASIS (2026-09-01, second iteration of this design): the
       first version used `bootstrap_rec3D`, which derives sigmas from
       the startup half maps. That exists for the case where no
       box-compatible sigmas CAN exist (the final rec, crop box ->
       native box); at the refinement box calc_pspec can estimate them
       directly from the particles, and its estimate conditions the
       euclid system far better (bgal startup base residual 0.23 with
       the half-map sigmas, 0.08 with calc_pspec's). Worse, refine3D
       then re-derived its own regardless, because
       `sigma2_stage_needs_bootstrap` was positional
       (`startit <= 1`) and never checked whether usable sigmas
       existed -- so the startup was regularized against sigmas the
       refinement discarded. That predicate now also requires the
       absence of a sigma star for the starting iteration, so one basis
       serves the whole workflow. Fresh runs are unaffected (no star
       exists, so the bootstrap still runs).
       The startup's `PCG_NU_STATS_FILE` is deleted afterwards: the NU
       replay controllers compare consecutive REFINEMENT iterations,
       and leaving the startup readout behind made the first
       auto-target comparison startup-vs-iteration-1 -- not like for
       like, and it took a real setpoint step off it (bgal: 15 -> 9).
   (c) The matcher can no longer match unfiltered references silently:
       when the nonuniform references are missing it flags the
       fallback, applies the FSC-optimal filter to the raw half maps
       instead of nothing, and hard-errors if there is no FSC either.
       This is the structural guard -- it protects every workflow, not
       just this one, and would have made the original defect loud.

6. **PfCRT matching-reference collapse -- SUPERSEDED HYPOTHESIS
   (2026-09-01, retained as a record of what was falsified).**
   Reproduced twice in refine3D_auto+pcg: iteration 1 matches with the
   SPHERICAL reference mask (the envelope does not exist yet) and is
   healthy (orientation overlap 0.875); iteration 2 is the first to
   consume the NU-evidence-envelope-masked references and collapses
   (overlap 0.236, in-plane distance 17 deg, shift increment 1.4 px
   avg / 11.3 max, FSC 0.500 6.47 -> 7.39 A, B-factor -79 -> -142).
   Mechanism: the NU evidence envelope OMITS DETERGENT BY DESIGN
   (PfCRT occupancy 8.9% of the spherical support), so the matching
   reference is protein-only while the particle images contain the
   micelle. Under the Euclidean objective the unexplained micelle
   dominates the cost, the assignment randomizes, and the
   re-centering shows up directly as the shift-increment jump. This
   is a MODEL-DATA MISMATCH, not overfitting and not a prior
   miscalibration -- soluble specimens (streptavidin, bgal) are
   immune because their envelope covers the whole particle.
   FALSIFIED by the automsk=no control, which reproduced the collapse
   identically (6.75/3.98 -> 7.76/4.93 with no envelope masking
   anywhere, versus 6.47/3.98 -> 7.39/4.93 with it). The reference
   mask was never the cause; see item 7 for the actual one. The
   matcher mask change this motivated was reverted. Retained because
   the reasoning about reference/image content consistency is still
   the right frame -- it simply pointed at the wrong mask.
   - (reverted) `prepare_matching_reference_mask` using the density
     envelope for matching references.
   - refine3D_auto gains a STARTUP BOOTSTRAP (user-directed):
     bootstrap_rec3D runs before any matching (unregularized pass ->
     sigmas from the half maps -> regularized pass with the
     refinement's own filtering settings), so the masks, the NU
     evidence products and the _nu_filt matching references all exist
     at iteration 1. Iteration 1 then matches the SAME kind of
     reference as iteration 10; previously the reference convention
     changed underneath an already converged alignment, which is what
     made the mismatch above fire as a cliff rather than a gradient.
     abinitio3D never needed this because its stage ladder introduces
     ml_reg / NU filtering / envfsc / automsk gradually, while the
     alignment is still loose.

8. **PfCRT refine3D_auto ROOT CAUSE -- COARSEST-BANK FLATTENING IN THE
   MATCHING REFERENCES (2026-09-02, confirmed by elimination; fix
   implemented).** With the startup bootstrap in place (1d0071e77) the
   collapse moved from iteration 2 to iteration 1 and tracked the first
   consumption of NU-filtered matching references in every topology.
   Control chain that isolated it:
   - automsk=no reproduced the collapse (mask exonerated -- and the
     earlier automsk=no falsification of the envelope hypothesis was
     confounded by the then-unfixed unfiltered-refs defect of item 7;
     both hypotheses pointed at reference/image content mismatch, both
     via the wrong mechanism);
   - rec_backend=gridding reproduced it (Q_NU replay, auto-lambda,
     evidence reuse, PCG solve support all exonerated);
   - `filt_mode=none automsk=no autoscale=no` was HEALTHY at iteration 1
     (overlap 0.981, angular moves 0.98 deg avg, shifts 0.17 px, FSC=0.143
     holding at 3.98 A, cFAR 0.73, B-factor -70) -- exonerating the
     euclid sigma cold start, the native box-300 geometry, nspace 20000,
     and prob_neigh, since only the reference recipe differed.
   Mechanism: on a detergent-solubilized specimen the NU bank puts
   72-73% of the spherical support in the coarsest 19.4 A bank (solvent
   clamp ~49% plus in-envelope null ~24%, i.e. ~45% of the density
   envelope interior -- the micelle). The matching reference then omits
   density every particle image contains; under the euclid objective the
   unexplained micelle dominates the cost, pose discrimination is lost,
   and the alignment scatters (3.93 -> 4.93 -> 5.86-5.97 A plateau read
   as convergence). abinitio3D is immune (nonuniform_lpset merged refs
   at box 160, where null ~= clamp only); compact soluble specimens are
   immune (nearly nothing inside their envelope goes null). Score-spread
   flatness (0.380 +/- 0.002) appears in healthy runs too and was never
   evidence.
   **First fix attempt -- FSC-OPTIMAL FALLBACK IN THE MATCHING PRODUCT
   (6e5135689) -- TESTED AND RETIRED (2026-09-02).** Seeding the
   `_nu_filt` matching references with the FSC-optimal filtered base
   pair wherever no finer bank label was positively selected did NOT
   rescue the alignment. The test ran at automsk=yes defaults, where
   the matcher still multiplied the references with the NU evidence
   envelope before reprojection -- the restored fallback content
   outside the ~10-19% envelope was erased again by the
   multiplication, which stood in every collapsing configuration.
   The multiplication, not the filter flattening, is the standing
   suspect; the fallback machinery was reverted.
   **Superseding policy (2026-09-02, user-directed).** Principle: a
   reference must never HARD-REMOVE density that is present in the
   particle images; density excluded by an envelope is DOWN-WEIGHTED
   through heavy low-pass filtering of the background (cisTEM
   precedent: envelope + heavy background low-pass does not break
   alignment). Three parts, both backends:
   - Reference-envelope multiplication REMOVED. The matcher's
     `mask_matching_reference` applies the spherical soft mask ONLY;
     `prepare_matching_reference_mask` and the evidence/density
     envelope multiplication (`zero_env_background` + `mul`) are
     deleted. No reference is ever multiplied with an envelope before
     reprojection.
   - automsk=yes now means: the filter-field BACKGROUND is the
     complement of the NU evidence envelope, derived from the SAME
     evidence pass (the setup unaries; no automask3D double-compute).
     `write_nu_evidence_envmask` gains `l_arm_background`: one call
     derives the envelope, writes the artifact, and arms the
     background clamp (`set_nu_solvent_envelope` with
     source='nu_evidence_envelope'; provenance -- and thus the frozen
     evidence identity -- records the clamp source). Background voxels
     take the coarsest bank candidate in the filter field, which both
     generates the matching references (heavy background low-pass,
     detergent down-weighted but present) and derives the Q_NU
     precisions (maximum lack-of-evidence outside the envelope) -- one
     field, both consumers. Armed before `optimize_nu_cutoff_finds` at
     all three sites (`build_nu_replay_evidence`,
     `filter_pcg_nonuniform_maps`, volassemble's
     `setup_nonuniform_filter`), so the existing pre/post-Potts clamp
     mechanics apply unchanged; the later plan-gated envmask writes are
     skipped in this mode (already written in-pass).
   - automsk=no keeps the conservative density-envelope solvent
     constraint of dev item 4 (on-the-fly automask3D), bit-identical.
   - The PCG SOLVE support is deliberately NOT the evidence mask: it
     stays the conservative density envelope from the lag-one
     reference (`build_pcg_state_support`), per the solve-support
     policy -- the evidence mask is too tight to constrain the solve.
   - Validation plan: PfCRT refine3D_auto at full defaults -- expect
     iteration 1 to hold ~3.9 A with high overlap now that the
     references explain the micelle at 19.4 A instead of omitting it;
     streptavidin rerun -- expect no regression (its evidence and
     density envelopes nearly coincide, and its background is
     essentially solvent).

9. **PCG+NU code-review response (2026-09-02,
   `pcg_nonuniform_code_review.md`, corrected implementation).**
   - P1 post-hoc filter: `filter_pcg_nonuniform_maps` no longer runs
     `nu_filter_vols` for ANY mode -- it is the evidence-derived scalar
     matching-lp handoff for both `nonuniform` and `nonuniform_lpset`,
     nothing else. No `_nu_filt`/`_nu_locres` products exist on the pcg
     backend; plain nonuniform matches the PRIMARY per-half Q_NU maps
     (matcher `l_pcg_qnu` path: primary maps authoritative, no fallback
     flag, no FSC-optimal or other extra filter). The evidence-phase
     setup-retention machinery in the strategy was removed as its
     consumer no longer exists (the nu_filter module API remains).
     NOTE the standing tension with the 2026-08-31 "matching-reference
     regression -- plain nonuniform topology" record, which retired
     exactly this raw-Q_NU-refs matching for gold-standard noise
     overfitting: circumstances have changed (evidence-envelope
     background suppression in-solve, auto-lambda/auto-target live),
     but the overfitting risk of independent per-half matching against
     each half's own in-band noise is not structurally removed --
     watch mid-resolution half-map divergence with a matched-noise
     0.143 tail on validation runs.
   - P1 evidence background: the solvent/background clamp in
     `build_nu_evidence_state` now encodes the COARSEST bank candidate
     (supports [1,0,0,...] -> band weights [0,1,1,...]) instead of the
     post-Potts null label whose zero supports suppressed all non-DC
     background detail. The constraint is installed before the Potts
     sweeps, fixed voxels participate as boundary conditions, and the
     readiness `null_fraction` is retained from the unconstrained broad
     spherical evidence field, counted only over that packed support, so
     neither the full-box exterior nor envelope size can invalidate its
     calibration diagnostic. The coarsest band's support is categorical and
     is not inferred by comparing a grid-quantized cutoff with its nominal
     Angstrom boundary.
   - P1 solve support: `build_pcg_state_support` is independent of
     `automsk`; the shared route now builds and installs the per-state
     density support exactly like the distributed route (base solves
     under envfsc=yes, regularized replay always). A valid start volume
     is a density source. If no reconstruction exists yet, the base solve
     necessarily bootstraps on the sphere for either `envfsc` value; its
     current merged pair then seeds the density-supported replay.
   - P1 FSC preproc: `evaluate_halfmap_pair` takes
     `l_pair_support_constrained` and skips the envelope+
     phase-randomization preprocessing only when the pair was actually
     density-constrained in the estimator -- never inferred from the
     backend name. PCG spherical-fallback pairs under envfsc=yes now
     get the ordinary envelope-corrected FSC.
   - P2 Q_NU mandatory: `validate_nu_replay_request` hard-errors any
     pcg+NU configuration without the active euclid ML replay and
     pcg_nu_lambda_rel > 0; strength-zero/P_tau remain valid only with
     filt_mode=none (bootstrap_rec3D pass 1 forces filt_mode=none and
     is unaffected). `pcg_mskfile` is likewise rejected on NU routes
     (development escape hatch isolated to filt_mode=none).
   - P2 refine3D_auto envfsc: now a guarded (genuinely overridable)
     default.
   - P2 _pproc: PCG skips all post-hoc mask multiplication in
     `postprocess_volume_from_files`, including derived `_pproc`/`_mirr`
     products. Non-PCG behavior is unchanged.
   - P3: `nu_evidence_envelope_masking.md` carries a SUPERSEDED banner;
     automasking/nonuniform policies updated 2026-09-02.
   - Runtime follow-up: a density-supported iteration-1 base pair exposed a
     numerical mismatch between the narrower solve support and broader NU
     evidence sphere. Exact zero/zero boundary voxels had entered the radial
     MAD whitening fit as zero-noise samples, and the explicit zero predictor
     could then overflow even though ordinary cross-half candidate unaries
     remained finite. `nu_objective_noise_profile` now excludes only those
     unobserved exact-zero pairs from scale estimation while retaining the
     spherical evidence domain; `nu_objective` computes standardized Huber
     losses in double precision with a high finite cap. The focused envelope
     test includes this hard-supported-pair geometry. Rebuild and runtime
     rerun are pending.
   - `nu_refine` follow-up: `refine3D_auto` keeps its default `yes`, while
     staged abinitio3D's heavily used `nu_refine=no` PCG state is preserved
     exactly (eight signal candidates, integer coordinates, unit candidate
     masses, four bands, raw-finest matching handoff). Adaptive PCG discovery
     alone is capped at FSC=0.143 plus two Fourier shells and uses the final
     ordered-label tie tolerance. Once extension candidates exist, their
     coordinates continue from the unchanged static ladder in Fourier-shell
     distance normalized by the finest static interval; normalized Voronoi
     widths remove candidate-count bias from soft band support. Adaptive
     matching uses the 5% supported cutoff plus the same two-shell FSC
     headroom. The `automsk=yes` evidence envelope is explicitly the
     preliminary static-bank boundary, fixed before adaptive challenges.
   - Observed-domain calibration (review follow-up, 2026-09-02 pm): the
     whitening fix above left every OTHER evidence statistic on the full
     sphere. With the base pair density-constrained (refine3D defaults,
     iteration 2 onward) about half of the spherical support is an exact
     zero/zero cluster at gap ~0, which pinned the lower-quartile null
     center, diluted the spatial beta, could collapse the temperature
     median into its fallback, and inflated the readiness `null_fraction`
     by support geometry (a small particle in a generous sphere could
     cross the 0.90 ceiling and hard-error). `setup_nu_dmats` now records
     a packed observation mask (`nu_observed_mask`, same exact-zero test
     as the profile); `build_nu_evidence_state` evaluates the null-bias
     center, beta, temperature, null/uncertain/band-support fractions on
     observed voxels only, freezes unobserved voxels at the explicit null
     in both Potts sweeps, and ships them with zero band support and unit
     uncertainty. The summary carries `observed_fraction` and the
     provenance string records `statistics_domain=observed_support`. The
     spherical NU support itself is unchanged (skill invariant).
   - Solve-support provenance: every shipped state volume now has a
     `<vol>_pcg_support.txt` sidecar (`solve_support=density|sphere`).
     The trailing bootstrap reads it for the lag-one FSC/evidence pair, so
     `evaluate_halfmap_pair` skips the envelope+phase-randomization
     preprocessing exactly when that pair was density-constrained in the
     estimator; a pair without a sidecar (imported) stays unconstrained.
     A bootstrap blend is constrained only if both contributions were.
   - Shell-walk flag renamed `l_require_margin` (was `l_tie_tolerant`):
     it demands a `nu_label_smooth_is_better` margin win, i.e. it is
     STRICTER than the raw `<` test, which the old name inverted.
   - Diagnostics `half_pair_parallel=` now reports whether the pair
     actually ran as concurrent sections (one empty half runs serially).
   - Shared route: a missing evidenced matching low-pass is logged and the
     previous `lp` rides, matching the distributed route.
   - `bootstrap_rec3D` imports `NU_AUTOTARGET_MIN/MAX` instead of
     duplicating the bounds.

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
