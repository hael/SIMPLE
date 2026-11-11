# Continuous local Cartesian refinement — done

Updated 2026-11-11. This records delivered mechanisms and observed validation,
not scientific acceptance of a width policy or resolution of the new review.
The **single active work list** is
[the todo record](../planned/refine3d_local_cartesian_priors_todo.md), containing
the remaining technical concerns and their response/acceptance plan.

Split at the maintainer's request from `refine3d_local_cartesian_priors.md`.
The original feature is `f68dab7b0` on parent `4965b3fe8`; the architectural
refactor was reviewed in the working tree against `13ad4833e`.
No commit or new runtime validation is implied by this documentation split.

## Delivered implementation

- Single-state Euclidean Cartesian joint refinement accepts the opt-in
  `cont_prior=none|shift|pose` controls. Shift penalties use fixed or all-active
  stored-increment widths; pose mode adds an explicit-width Matrix-Fisher
  rotation penalty. The current penalty re-centres on each transaction's seed.
- Numerical terms, native/cropped-pixel scaling, particle-power access,
  set/clear lifecycle, bounds, rollback and parameter/worker propagation were
  implemented with existing numerical and workflow regression coverage.
- Optimization uses the penalized objective; committed scores and sigma2
  residuals retain their data-only meaning. Prior potency, convergence effects
  and statistical calibration
  remain open in C2-C8 of the todo record.

## Fixed shift-prior command-line units

The maintainer requires the same command-line convention as polar refinement:
positive `cont_prior_shift_sigma` values are per-axis standard deviations in
**Angstroms**. The parameter remains in Angstroms through construction and
worker serialization. Pass preparation converts once to native-pixel variance,
`(cont_prior_shift_sigma/smpd)^2`, using native `smpd`, independently of cropping.
The optimizer retains its existing native-to-cropped-pixel conversion.

Zero still selects the all-active automatic estimator and its native-pixel
floor; prior-off behavior is unchanged. Historical native-pixel CLI values in
the archive below describe the old snapshots. To reproduce an old fixed width,
use `sigma_A = old_sigma_native_pixels * smpd`; for example, 2 pixels at
1.5 A/pixel is now entered as 3 A. Historical measurements are not rewritten.

Unit-conversion regression coverage exercises 3 A at native samplings 1.5, 2
and 3 A/pixel, at full and cropped boxes, retaining the raw command-line value
across repeated preparation. It also covers automatic-floor preservation,
fixed widths below that floor, and clearing the disabled prior. Build/runtime
validation of this unit correction is recorded separately from the earlier
refactoring evidence below.

Observed for this correction on macOS, GNU Fortran 16, the existing Release build:

- `make -C build -j8` and complete UI JSON validation passed without compiler
  warnings/errors; `make -C build install` passed after the gate.
- Focused `unit_cart_align3D suite=pose_strategy` passed 100 assertions,
  including the new Angstrom conversion cases. All six static checks and the
  full fast gate passed, 15/15 in 58.23 s.
- `lib_cart_align3D` passed in 23.80 s. The existing
  `cont_refine3D_1jxy suite=nu_handoff` passed in 19.3 s against the refreshed
  installation. This was the focused integration case, not another full
  continuous-workflow run or a real-data prior-quality comparison.
- Navigation was regenerated; changed doc links and whitespace were checked.
  Logs and the pre-change snapshot are in `/private/tmp/simple-prior-angstroms/`.
  No new Linux/BOX validation, commit or push.

## Delivered architectural refactoring

- Cartesian names describe pose updates and a continuous followup pass:
  `get_pose_update`, `l_cont_followup`, `calc_cont_pose_stats`, and corresponding
  pass objects/callers. Benchmark `suite=local` replaces `suite=polish` in code
  and help; the public `pose_cont` control and stored project fields retain
  their contracts.
- Four optimizer kernels became private module procedures in the production
  submodule `src/main/cftc/simple_cartft_pose_opt_system.f90`. Its sibling tester
  retains access. The shared mean-square shift helper is private, with its
  closed-form behavior covered through the production width API.
- Declaration/assignment layout was adjusted and navigation regenerated.
  The attached review identifies residual alignment errors; those remain C9,
  rather than being claimed complete here.
- The numerical source-equivalence and build/test evidence below were observed.
  No adaptive cone, posterior exporter or posterior comparison harness was
  implemented in that refactor.

## Latest observed validation before this split

Observed locally on macOS, GNU Fortran 16, the existing Debug build:

- `make -C build -j8` and complete UI JSON validation passed. All six static
  gate checks passed; final fast gate 15/15 passed in 145.07 s after the last
  visibility cleanup (the earlier run also passed, 145.75 s). The 60 s overrun is
  advisory; sandbox CPU-count detection selected serial suite execution.
- `make -C build install` passed after the final gate. The new production numerical submodule is
  intentionally untracked until the maintainer includes it with these changes.
- Parent, numerical submodule and sibling tester compiled/linked at `-O3`;
  the existing `pose_optimizer` suite passed 118 assertions against those
  optimized objects. External `USE` imports of each of the four optimizer
  helpers and the mean-square shift helper were rejected by the compiler.
- All eight moved kernel bodies match `13ad4833e` after whitespace/module-prefix
  normalization. All eleven changed production logic files preserve tokens after
  the documented API renames/extraction/accessibility changes, excluding messages
  and whitespace (`/private/tmp/simple-cont-review/source-equivalence.json`).
  This is a source-equivalence check, not a full Release workflow identity run.
  The additional probability-utility edit changes only the mean-square
  helper's accessibility; its implementation is byte-identical to the baseline.
- `lib_cart_align3D` passed in 91.08 s; `cont_refine3D_1jxy` passed in
  1462.63 s (135 metric/check rows, no failed gate rows). These exercised the
  numerical extraction and terminology changes, including shared/distributed
  agreement, local continuation, followup sampling/reconstruction, in-place
  refinement and NU handoff. The later visibility-only mean-square helper
  cleanup and removal of its duplicate direct assertion were rebuilt and are
  covered by the final fast-gate rerun; the long workflow was not repeated.
- Benchmark `suite=list` reports `local, capture`. Registry, FLEX DAG, source
  description, documentation links and whitespace checks passed. Repository
  navigation was regenerated as required by `.codex/AGENTS.md`.

Logs and compiler/source probes are retained in `/private/tmp/simple-cont-review/`.

No new Linux/BOX or real-data quality/speed qualification. Historical Linux
evidence in archived section 3 does not qualify these working-tree changes. Adaptive
restraints and the proposed posterior comparison harness remain unimplemented.
No commit or push.

## Documentation split verification

On 2026-11-11, the combined record was replaced by this `_done` archive and
the `_todo` plan at the maintainer's explicit request. The original sections
1-3 are retained below, unchanged except for heading depth. All nine concern
headings and six posterior follow-up topics were checked against the supplied
review. Changed local links/fragments, fences and whitespace were checked;
the index, local algorithm and local policy now point to the two records.
All 865 Fortran source files matched their pre-split hashes. This was a
documentation-only change; the runtime results above were not rerun.

## Archived contract, implementation plan and delivery history

The following preserves the original sections 1-3, with headings nested one
level deeper. It is a historical record, **not a second active work plan**.
Future-tense work packages, feature-only commands, superseding progress entries
and old approval wording describe their original snapshots. All still-open
work is tracked in the todo record.

In particular, the historical curvature comparison and the
“identifiability” promise are under review in C4 and C5. Historical passing
checks do not establish restraint potency, safe natural convergence or calibrated
posterior uncertainty. Do not read an old FINAL/COMPLETE statement as closure
of the newly attached concerns.

### 1. Specification: what is proposed

#### Purpose and scope

Add opt-in update regularization to single-state `refine3D_local` with `objfun=euclid`:
shift prior first, then Matrix-Fisher rotation prior after separate review. Preserve
the Cartesian Levenberg-Marquardt (LM) stage configuration and coupled solve of all
five parameters (three rotations, two shifts).

Prior support is limited to the production joint route. Reject
`cont_route=shift_then_joint` with `cont_prior /= none`; do not extend `refine_shift`.
Retirement of that internal test route is separate work. Prior-off behavior remains
covered on existing routes.

Assume late-stage refinement with a useful stored pose. Regularization around each
transaction's seed can suppress noise-driven motion but impede correction of wrong
seeds. It neither escapes local minima nor limits drift from the original poses.

This living note replaces and corrects `cartesian_posterior_bites.md`; older
neighborhood/cycle plans do not govern it. Polar search, top-k reconstruction,
active-state search, posterior sampling and Laplace output remain out of scope.

#### Public control

Use the maintainer-accepted `cont_prior` control below. Hans's review permits shift
work once its required design revisions are incorporated; scientific benefit remains
to be established by tests.

| Value | Meaning and availability |
| --- | --- |
| `none` | No prior; default. |
| `shift` | Shift penalty; first implementation. |
| `pose` | Shift and Matrix-Fisher rotation penalties; requires explicit `cont_prior_rot_sigma` in degrees. |

`cont_prior_shift_sigma=0` uses the all-active estimator; a positive value fixes the
per-axis shift standard deviation in native pixels. `cont_prior_rot_sigma` is a
positive, explicitly supplied per-axis standard deviation in degrees for `pose`.
Both controls were approved by the maintainer during implementation; no rotation
width is inferred from the LM bound or chosen as a production default.

This selects regularization, not active LM parameters. Advertise it in
`refine3D_local`, and permit explicit use in compatible base `refine3D` calls:
single-state Euclidean continuous refinement, not auto polishing. The maintainer
accepted this scope clarification. Default-off behavior prevents automatic activation
in other callers.

| Alternative | Benefit | Cost or risk | Decision |
| --- | --- | --- | --- |
| Reuse polar `shift_prior` | No new key. | Defaults to `yes`; changes its meaning and can activate unrelated Cartesian callers. | Reject. |
| Commander-only derived flag | No public key. | Hidden behavior; a local flag does not survive child parameter construction by itself. | Reject as the sole control. |
| Explicit `cont_prior` | Clear opt-in behavior; transferable through existing parameters and command lines. | Requires UI/parser and propagation coverage. | Selected, default off. |

#### Mathematical sources and implementation provenance

These references document mathematical foundations; they are not retrospective
claims that external implementations were copied or consulted during development.

- **Gaussian shift penalty:** the negative log of an isotropic Gaussian is
  `||delta_s||^2/(2*sigma_s^2)` up to a constant; its gradient and Hessian follow
  by differentiation. Background: C. M. Bishop, *Pattern Recognition and Machine
  Learning* (2006), [author/publisher resource](https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/).
- **Rotation model:** T. Lee, *Bayesian Attitude Estimation with the Matrix Fisher
  Distribution on SO(3)*, IEEE TAC 63(10), 3377-3392 (2018),
  [doi:10.1109/TAC.2018.2797162](https://doi.org/10.1109/TAC.2018.2797162);
  [author manuscript, section II-A, equation (5)](https://arxiv.org/html/1710.03746v3).
  Substituting `F=kappa*R_seed` into `p(R) proportional to exp(trace(F^T R))`
  gives the isotropic trace penalty. Its matrix-difference form, small-angle width
  matching, right-increment derivatives and Gauss-Newton curvature are derived
  for this implementation; the paper does not supply our crop/power scaling.
- **SIMPLE reuse:** `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`:
  existing `shift_prior_sigma2` supplies the mean-square formula and floor.
  The extracted helper preserves polar population and rounding. Existing
  `src/main/cftc/simple_cartft_pose_opt.f90` supplies joint LM, coordinate scaling,
  Cholesky, bounds and transaction behavior. The all-active Cartesian population
  is our experimental policy, not a literature-validated estimator.
- **Local derivation:** `T_cart=P/4`, crop conversion, penalty weights and the DC
  approximation follow SIMPLE's Fourier counting/noise and normalized-loss
  conventions (Numerical contract below), not a borrowed universal formula.
  Independent scaling/derivative tests support the implementation; they do not
  establish scientific reliability or justify the quarter-wavelength floor.

#### Numerical contract

Let `L` be normalized data loss, `P` whitened particle power, `c = box_crop/box`, and
`delta_s = s - s_seed` in cropped pixels. Production observations already include
the stored shift, so their seed increment is zero. Lower-level seeds may be nonzero.
Always subtract the actual seed and freeze it through all stages.

Here sigma2 is the variance per real/imaginary Fourier-noise component. Each conjugate
pair is one independent complex observation counted twice in the full-plane residual
sum. Its negative log-likelihood is one quarter of that sum, giving `T_cart = P/4`.
Independently test this convention, normalization and active-band/DC handling.

Accepted DC exception: when the active band includes the real zero-frequency pixel,
`P/4` gives its likelihood contribution half the required weight. Retain this
approximation and the existing data objective for the first experiment. Test
conjugate-pair scaling and the DC exception separately; state whether DC is active.

The Gaussian shift penalty in loss units is:

```text
F_s = w_s * dot_product(delta_s, delta_s)
w_s = 1 / (2 * sigma_s^2 * c^2 * T_cart) = 2 / (sigma_s^2 * c^2 * P)
g_s = 2 * w_s * delta_s
H_s = 2 * w_s * I_2
```

Here `sigma_s` is a per-axis width in native pixels. Do not apply `prob_beta`:
`T_cart` converts units and is not a sampling schedule.

| Requirement | Contract |
| --- | --- |
| Data separation | Calculator loss, score and shell residual remain data-only. Keep separate data and penalized objectives in the optimizer. |
| Acceptance | Compare penalized objectives. The seed penalty is zero, so final acceptance requires data improvement greater than the final penalty. An intermediate step may increase data loss while lowering the total. |
| Identifiability | Check the data-only Hessian before adding prior curvature. Split the check from LM system construction, reusing existing coordinate scaling, active-parameter masking and Cholesky logic. Reject an unidentifiable data system; the prior cannot rescue it. |
| Curvature and step control | Add the exact prior Hessian after data `curvature_scale`, without multiplying the prior by that heuristic. Stationarity uses the combined gradient; predicted reduction uses combined derivatives; actual reduction uses the penalized objective. Gain ratio, damping adjustment and predicted-reduction stopping use these combined quantities. |
| Bounds and lifecycle | Preserve `trs`, `athres_cont`, stage configuration, rollback, damping algorithm and tolerance constants. With an enabled prior, feed penalized quantities to optimizer step/stopping checks; keep workflow motion-based convergence unchanged. Apply the penalty once at every evaluation, including cached terms. |
| Disabled path | `cont_prior=none` preserves existing results bit for bit on the same build and fixture, including statuses, scores and motion fields. |
| Objective | Enable the prior only for Euclidean loss; do not silently apply its scale to CC. |
| Precision | Use `real(dp)` for all new widths, reductions, powers, crop factors, weights, penalties and derivatives. Promote operands before arithmetic. Preserve existing double-precision LM accumulation and solving. |
| Ownership | Numerical terms belong in `cftc`; strategies configure them. No new project I/O in workers, duplicate parameter objects, new CMake targets or compiler-wide precision flags. |

#### Published quantities

| Quantity | Meaning with an enabled prior |
| --- | --- |
| `get_objectives`: `objective_before/after` | Data-only at transaction input/returned pose; preserve existing score callers. |
| `corr_cart`, shell residuals and sigma contributions | Data-only. |
| `get_stage` objective values | Penalized optimization objectives; document and label this meaning in consumers. |
| Existing benchmark objective fields | Retain data-only meaning; add explicitly named penalized fields only where consumed. |
| `pose_cont_improved`, `pose_cont_status` | Penalized transaction acceptance/status. |
| Motion/convergence summaries | Actual pose movement and resulting status; retain existing workflow rule. |

The transaction seed penalty is zero; at the returned pose the difference between
penalized and data-only loss is its penalty. Audit `get_stage`, `get_objectives` and
work-record consumers together; update labels and assertions that previously equated
stage and transaction objectives. Preserve prior-off identity. This reporting policy
is accepted; exact benchmark columns follow their existing consumers, without adding
unused getters or fields.

#### Experimental shift-width policies

Policy in brief: the shift width is a soft leash around each LM run's starting
pose. A smaller width gives a stronger penalty; it is not a hard motion limit.
Automatic width uses all active records' latest displacement magnitudes, interpreted
as per-axis variance under zero-mean, isotropic 2D motion, with a one-native-pixel
or quarter-matching-wavelength minimum. This limits tightening as updates shrink.
The width stays fixed within each run and is re-estimated at the next pass; a
positive explicit shift sigma instead selects a fixed width. These assumptions
and the inherited minimum remain experimental.

The first adaptive candidate uses **all active particles' latest stored shift
increments**, not the current worker's partition. Production policy will be selected
from comparative evidence. For `N` active records its shared mean-square formula is:

```text
sigma_floor = max(1 native pixel, lp / (4 * smpd))
sigma_s^2   = max(sigma_floor^2, sum(shincarg_i^2) / (2*N))
```

The factor two converts mean squared 2D displacement to per-axis variance. Use the
matching low-pass, native sampling and double-precision reduction. Freeze the width
before updates. Fractional updates leave increments from different iterations: this
is not a strict previous-iteration estimator. Include finite, nonnegative values,
including zeros. The record initializes `shincarg` to zero, so absent history and
measured zero motion cannot be distinguished. All-zero history uses the floor;
report floor selection without claiming measured stability. A negative, NaN or
infinite increment in an active record stops with the particle identified. Do not
silently discard invalid records or change the population; inactive records do not contribute.

Use existing fatal-error and queue-failure handling. Width validation runs serially
before particle matching: report particle ID and stop nonzero. The queue controller
already detects failed/lost parts and stops the master (`src/utils/qsys/simple_qsys_ctrl.f90`:
`stop_on_failed_parts`); reuse it without a new failure protocol. Existing
`simple_qsys_ctrl_tester.f90` covers failed local parts/statuses. After implementation,
add/run a focused Linux distributed case asserting that an invalid active increment
causes master failure before merge, without waiting for successful completion or
accepting partial output. This is future coverage, not a test available to run now.

The scientific reason for this candidate population is to estimate typical latest
motion across the active project, rather than let a changing sample define the width.
This trades freshness for population coverage: unsampled particles retain old values,
and a polar-to-Cartesian handoff can mix discrete and LM increments. Cartesian and
polar priors may therefore have different widths on the same field. This rationale
and its scientific benefit need comparative evidence. Preserve polar
semantics. Scalar reset/preparation is specified in the ownership table.

Use fixed-width fixtures to isolate penalty mathematics. Adaptive widths may shrink
as the prior suppresses motion; test tightening and premature convergence. The floor
alone does not establish safety.

Agreed experimental comparison: prior off, adaptive width, and fixed standard
deviations of 1x, 2x and 4x the initial width floor. Freeze each fixed width throughout
its run; match starting poses, references, random seeds and iteration budget. These
are experimental settings, not production defaults. Derive numerical tolerances from
analytic formulas and finite-difference accuracy, not candidate results.

Production widths must follow an explicit statistical model and a sampling/resolution
rationale, evaluated across datasets. The beta-gal comparison checks behavior; it
does not determine a dataset-specific production constant. The mean-square formula
has a per-axis displacement interpretation; the one-pixel/quarter-wavelength floor
is an inherited safeguard whose production justification remains open.

A completed-sample median-based width is an additional proposed comparator, not yet
a specified estimator. Before testing it, define its conversion to per-axis variance
and capture the completed sample before next-sample markers overwrite its identity.
Do not substitute median displacement directly into the variance formula. Median
limits minority outliers but can hide difficult tails and cannot protect against an
iteration where most particles move incorrectly. Compare outlier sensitivity, stale
history and bad-iteration behavior; do not choose a production winner before testing.
This comparator's definition is a later experiment-design gate, not a blocker for
the first shift-penalty and fixed/all-active tests.

#### Rotation choice and width policy

Use an isotropic Matrix-Fisher penalty in the later rotation step:

```text
F_r = (w_r/2) * ||R - R_seed||_F^2
w_r = 1 / (2 * sigma_r^2 * T_cart)
```

For rotation matrices this equals `w_r * (3 - trace(R_seed^T R))`, approaching
`w_r * theta^2` at small angles. Matrix differences avoid trace cancellation and
give zero at rounded seeds. Use consistent right-increment derivatives and residual
Gauss-Newton curvature, not an assumed exact Hessian away from the seed.

| Alternative | Benefit | Cost or risk | Decision |
| --- | --- | --- | --- |
| Squared geodesic angle | Direct angular-distance interpretation. | Requires consistent rotation-log derivatives and curvature. | Not selected. |
| Squared angle with constant curvature | Simple local approximation. | Constant `2*w*I` is not its exact Hessian away from the seed. | Not selected. |
| Matrix-Fisher via matrix differences | Smooth residual formulation; no rotation logarithm; compatible with joint LM. | Width needs a stated convention; Gauss-Newton remains an approximation. | Selected. |

Start with fixed experimental widths: `sigma_r` is small-angle per-axis standard
deviation in radians, not RMS rotation magnitude. Compare recovery and stopping
against prior-off on known pose errors; values and margins remain open. `dist` and
`dist_inpl` are symmetry-dependent summaries, not rotation-vector components.
Defer adaptation until geometry, symmetry handling and an independent floor are
established. Do not reuse same-image posterior covariance as the next prior.
Rotation derivatives use physical radians. Apply `rotation_scale` through the existing
LM coordinate transformation after forming prior derivatives; do not encode that scale
again in `sigma_r` or the prior weight.

#### Acceptance and evidence

Use existing SIMPLE testers, failing assertions and fixed fixtures/seeds. Before
numerical implementation, freeze tolerances from independent formulas or pre-change
fixtures. Freeze scientific margins before acceptance evaluation, as specified below.
These are implementation requirements. See the progress record for completed checks and pending Linux validation.

| Coverage | Independent check and owner |
| --- | --- |
| Prior off | Exact baseline pose, status, score, residual and motion fields across existing optimizer stages; `simple_cartft_pose_opt_tester`. |
| Shift mathematics | Closed-form penalty, gradient and Hessian; finite differences at nonzero seeds; tight/loose widths; optimizer tester. |
| Scale and units | Known whitened power and Fourier-noise fixture; native/cropped pixel equivalence; calculator and optimizer testers. |
| Transaction behavior | Data gain smaller than penalty rejected; hard bounds, rollback, cached terms and joint-route seed anchors preserved; unsupported shift-then-joint/prior combination rejected; optimizer/parameter testers. |
| Identifiability | A fixture singular specifically in an active shift direction remains rejected even when its shift prior would make the system positive definite. Preserve prior-off results and active-parameter masking; optimizer tester. |
| LM model and reporting | Independently check unscaled prior curvature addition, combined predicted/actual gain and stopping inputs, plus data-only getters versus penalized stage records; optimizer tester and benchmark consumers. |
| Score and sigma | Same pose gives the same data-only score and residual with prior on/off; `simple_strategy3D_cont_tester`. |
| Width lifecycle | Hand-calculated variance; mixed zeros retained in the denominator; all-zero floor; negative/NaN/infinite active records rejected with particle ID; inactive records excluded; shared/distributed partition agreement; enabled/disabled reset before early return. Owner: `src/main/strategies/search/simple_strategy3D_cont_tester.f90`, existing `pose strategy` sub-suite of `unit_cart_align3D`. |
| CLI propagation | Default off; supported values and objective/scope validation; identical effective settings through shared and distributed child construction; existing parameter/UI tests. |
| Distributed failure | After implementation, an invalid active increment fails a worker and the master before merge; extend an existing distributed integration fixture, reusing queue-controller failure coverage. |
| Precision | Representative small/large finite power and width values; operand promotion before squaring; finite results and independent expected scaling; calculator/optimizer tests. |
| Later rotation | Matrix-Fisher residual and analytic gradient against finite differences; residual Jacobian and Gauss-Newton block checked separately; nonidentity seed. |
| Scientific benefit | Compare prior-off, fixed-width and adaptive-width cases: known-truth pose-error median/P90/tail, FSC and width evolution, then an aligned real project. Use matched iteration budgets and ordinary stopping to expose premature convergence. Unit correctness and the width floor alone do not establish scientific reliability. |

Report control, width and active count per pass. Reuse accepted/rejected and bound-hit
statistics; add penalty summaries through their owner. Smaller motion/loss alone
does not prove benefit. Select fixtures and practical run budgets before experiments.

Identify candidate scientific margins through exploratory tests on calibration
datasets. Declare calibration and held-out datasets (or independent realizations)
before tuning. Record the rationale, then freeze settings, margins and metrics before
independent validation. Do not count exploratory results as acceptance evidence under
thresholds chosen afterward. No separate approval of margins by Hans is required by
his written review; report the resulting evidence and limitations for maintainer review.

### 2. Implementation plan

#### Existing execution path and ownership

The local commander launches `refine3D`; its strategy prepares references and drives
shared/distributed iterations. The matcher reads the Cartesian band, selects particles,
prepares state/batches and calls `strategy3D_cont`. One optimizer transaction per
particle produces pose, data-only score and residual. Existing reconstruction and
sigma updates follow matching.

| Repository path / owner | Planned work |
| --- | --- |
| `src/main/params/simple_parameters.f90` | Declare `cont_prior` and its default. |
| `src/main/params/simple_parameters_parse.f90` | Bind the character option in `bind_input_carg`. |
| `src/main/params/simple_parameters_phases.f90` | Validate values and compatible objective/mode; derive effective settings. |
| `src/main/ui/simple/simple_ui_refine3D.f90` | Expose the approved local control with explicit choices and default. Generated argument/default files are regenerated, not hand-edited. |
| `src/main/commanders/simple/simple_commanders_refine3D.f90`: `exec_refine3d_workflow` | Preserve the control in the base-refinement child; keep orchestration and public-scope validation here. |
| `src/main/strategies/search/simple_strategy3D_alloc.f90`: `clean_strategy3D`, `prep_strategy3D` | Reset the separate Cartesian scalar; derive the experimental all-active width before the Cartesian return. Preserve polar allocation/estimation. |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`: `shift_prior_sigma2` | Reuse/extract the mean-square formula with caller-selected populations; preserve existing polar interfaces, population, floor semantics and rounding behavior, checked against baseline fixtures. |
| `src/main/strategies/search/simple_strategy3D_cont.f90`: `srch_cont` | Configure/clear the particle prior from pass state, crop scale and slot power; report data-only results. |
| `src/main/cftc/simple_cartft_calc.f90` | Add the slot-power getter near existing getters; leave normal terms and shell residuals unchanged. |
| `src/main/cftc/simple_cartft_pose_opt.f90`: `refine_pose`, `refine_joint` | Mirror set/clear lifecycle patterns. Freeze anchors; add the joint-route penalty after data evaluation and identifiability checks. Keep cached data terms unpenalized until used. Do not extend `refine_shift`; reject its test route with an enabled prior. |

Extend `run_all_strategy3D_cont_tests` in the named strategy tester for population
width and reset coverage. Exercise the real `prep_strategy3D`/`clean_strategy3D`
lifecycle on an in-memory fixture, with independently calculated widths; do not
duplicate the estimator in the fixture. Reuse its existing suite registration.

Follow the existing CLI chain: parameter declaration/default -> `add_char` binding
-> phase validation -> UI choices -> JSON validation. Initially declare/expose only
`none|shift|pose`; require explicit positive rotation width for `pose`. Existing CMake rules regenerate
`simple_args.f90`, `simple_varlist.txt` and `simple_ui_default_values.f90` when their
inputs change (`src/CMakeLists.txt:10-58`); no generator or CMake changes are expected.
Check `scripts/check_cli_chain.py` after implementation, plus Linux parameter/UI tests
for default `none`, accepted `shift`, invalid values and missing rotation width and JSON metadata.

Use one shared mathematical formula for shift-prior variance: the caller explicitly
selects the increment population; the shared calculation applies the mean-square
factor and floor. Keep Cartesian inputs, reduction and result in double precision.
Use a separate `real(dp)` scalar, `s3D%cont_shift_prior_sigma2`, for Cartesian width.
Reset it to zero in `clean_strategy3D` and calculate it per pass before the Cartesian
early return. Retain existing polar scalar types and optimizer interfaces; do not
promote the whole polar chain for this feature. Sharing the formula must preserve
polar behavior, verified by regression tests. This shift-prior variance is distinct
from image-noise sigma2.
Shared-memory parameters/project now persist across iterations
(`src/main/strategies/parallelization/simple_refine3D_strategy.f90`: `inmem_execute_iteration`).
Recompute width every pass; never rely on parameter reconstruction to clear it.
Verify equivalent loaded, all-active fields in both execution modes before updates;
do not use partition-only arrays, reread projects or add a scheduling layer.

Accepted identifiability mechanism: separate the data-only check from
`build_pose_lm_system` in `src/main/cftc/simple_cartft_pose_opt.f90`. At each joint
linearization, scale and mask the data Hessian using the existing operations, then
test it with the existing Cholesky logic. On failure retain the existing rejection;
on success add prior derivatives and construct the damped solve. Reuse the scaling
and factorization code without duplicating the check's implementation. Keep cached
data derivatives separate from penalized derivatives. Prior-off must preserve the
existing result. Stationarity follows the accepted step-control policy below.

Accepted step construction: `H_total = H_data_scaled + H_prior`, where
`H_data_scaled` is the calculator output after `curvature_scale`. Combine data and
prior gradients/objectives consistently, then apply existing LM coordinate scaling,
active-parameter masking and damping. Stationarity uses the combined scaled gradient;
predicted reduction uses combined derivatives, actual reduction uses penalized loss,
and their ratio drives damping. Apply `LM_PREDICTED_RELATIVE_TOL` to the penalized
objective scale. Keep existing constants and algorithms; the preceding data-only
identifiability check and published data-only scores/residuals remain separate.

Carry `cont_prior` on the existing command line: the local commander changes `prg`
to `refine3D`, distributed setup serializes it with `gen_job_descr`, and workers
parse their parameters normally. No origin flag or new worker transport is needed.
Validate compatible mode/objective/state count before launch and in worker parameter
construction. Strip the search-only control from non-search child commands using
the existing cleanup path in
`src/main/strategies/parallelization/simple_refine3D_strategy.f90`:
`strip_refine3D_search_only_args`; verify reconstruction/finalization callers use it.

#### Width-policy source evidence (rechecked 2026-11-09)

The builder loads the whole project (`src/main/simple_builder.f90`: project read),
but marker mutation matters:

| Path | Marker state before width preparation |
| --- | --- |
| Fractional Cartesian, distributed | Master selects globally and publishes the field in `distr_execute_iteration`; workers reproduce that sample. Markers agree, but selected particles still hold older increments before this search. |
| Full Cartesian, distributed | `master_draws_sample3D` is false; each worker selects its partition before `prep_strategy3D`. `sample4update_all`/`incr_sampled_updatecnt` advance only selected rows. Latest markers can therefore select different partitions in different worker copies. |
| Shared Cartesian | Selection also precedes preparation; the latest marker describes the newly selected set, whose new increments have not yet been computed. |
| Polar probabilistic | Master selects globally; workers reproduce it. The existing width is partition-independent on this path; this work does not change it. |

Source paths above: `src/main/strategies/parallelization/simple_refine3D_strategy.f90`,
`src/main/strategies/search/{simple_strategy3D_matcher,simple_matcher_smpl_and_lplims}.f90`,
and `src/main/ori/simple_oris_sampling.f90`. `get_all_sampled` selects markers above
the field maximum minus 0.5 (`src/main/ori/simple_oris_getters.f90`); it does not
establish that increments came from the last completed pass. Non-polish Cartesian
search writes `shincarg` only for processed particles (`simple_strategy3D_cont.f90`: `oris_assign_cont`).

Experiment candidate: all-active widths for population coverage and marker independence,
with the staleness/mixed-history limitation stated above. A strict last-completed-
sample alternative needs its population captured before selection, consistently on
the project owner, then passed to workers. Merely calling the existing getter at
the current preparation point does not provide that contract. Test both full and
fractional marker transitions; scientific advantage is not established by this trace.

#### Current implementation status

| WP | Source implementation and remaining work |
| --- | --- |
| 1 | Independent main log passed 80 assertions; its six printed baseline records match the earlier feature log exactly. Latest feature Cartesian regression passed 1,310 assertions before final cleanup. |
| 2 | Shift mathematics, lifecycle, LM integration and numerical fixtures implemented; Linux tests passed before final test deduplication. |
| 3 | Control, width/reset, particle setup and queue-boundary failure test implemented; Linux tests passed. This does not test an external scheduler. |
| 4 | Existing benchmark extended with six matched prior cases and median/P90/maximum pose-error metrics. Dataset runs, scientific margins and production settings remain evidence work. |
| 5 | Matrix-Fisher penalty integrated with explicit width, lifecycle and transaction tests; Linux tests and smoke benchmark passed. |

#### Linux validation and calibration commands

After the incremental build, run `simple_test_exec test=unit_cart_align3D` and
`scripts/run_fast_gate.sh build`. The `pose_strategy` suite includes the bounded
worker-to-master failure test; its expected child failures are asserted by the parent.
The test transfers a real worker exit status through the production queue boundary;
it does not submit an HPC job or exercise a particular external scheduler.

Feature-only experiment: the following command requires the feature worktree's
`suite=priors` benchmark extension; it is not included in this master delivery.
Run it in a new directory using that feature build:

```bash
set -o pipefail
SIMPLE_SEED=20261109 simple_test_exec test=cont_refine3D_bench suite=priors \
  objfun=euclid cont_prior_rot_sigma=5 maxits=1 nptcls=100 nthr=8 2>&1 | tee priors.log
```

`5` degrees is an explicit experimental example, not a production recommendation.
Use `nparts=2` to repeat the benchmark through the distributed matcher. The six
cases share observations, references and seeded perturbations; the smoke command uses one outer
refinement/reconstruction pass. Omit `maxits=1` for the six-iteration calibration
budget, with ordinary workflow and LM stopping. Retained `metrics.tsv`
reports median/P90/maximum errors, data score and accepted fraction. The smoke run does not establish multi-iteration scientific benefit or FSC acceptance.
Select further datasets and freeze margins before independent validation.

#### Work packages

1. **Prepare implementation:** re-check upstream source, capture baseline fixtures
   and numerical thresholds. Numerical support is the delivery goal; scientific tuning is optional follow-up.
2. **Implement shift mathematics:** power getter, prior lifecycle and shared helper;
   run optimizer/calculator coverage from the acceptance table first.
3. **Integrate pass policy:** control, population width and particle setup through
   the owners above; run strategy, lifecycle and command-propagation coverage.
4. **Evaluate shifts:** fixed-width cases, then population policy; compare matched
   budgets and ordinary stopping. Default changes require a separate evidence-based decision.
5. **Add rotation after review:** extend the helper with Matrix-Fisher and fixed-width
   tests; enable `pose` only when tested. Adaptive angular width remains deferred.

#### Precision audit: retained boundaries and follow-up

Power, objective, gradient, Hessian and LM solve/acceptance already use double precision.
This limits further arithmetic loss but cannot recover single-precision input detail.
The source audit established no observed overflow.

| Boundary in inspected source | Finding and disposition |
| --- | --- |
| `simple_cartft_calc:set_sigma2` (467), `fill_ptcl` (1167) | Inverse noise scaling, whitening products and particle/transfer storage are single precision. Extreme scaling can overflow before double accumulation. Retain initially; audit finite/range behavior in a separate measured change. |
| `simple_cartesian_fourier:gather_packed_bspline3_window_grad` (143) | Coordinates, weights, contractions and derivatives are single precision. Small-step sensitivity can be limited. A double-arithmetic gather could still read single-precision volumes; measure cost and accuracy separately. |
| `simple_cartft_calc:shift_phase_tables` (1273), optimizer bound conversion (151) | `real(PI,dp)` promotes an already rounded single-precision constant. `DPI` exists. Changing this affects the disabled baseline; keep separate from the prior patch. New prior arithmetic uses double-precision constants. |
| `simple_matcher_ptcl_batch:prep_cart_batch_of_build` (184), optimizer pose publication (217) | Stored shifts, crop preparation and published poses are single precision. Existing publication narrows before shift rescaling. Retain existing interfaces initially. |
| `simple_cartft_calc:sigma_from_shell_residual` (1124), strategy score publication (191) | Double-precision results are narrowed into existing sigma/score storage. Preserve their data-only meaning; audit conversion ranges separately. |

Retain finite/positive width and power checks and data-conditioning guards. Preserve
storage and default real kinds. Add negative-log penalties directly; existing `exp(-L)`
remains a data-only score, not posterior normalization.

#### Follow-up after the current delivery

Add prior support for `objfun=cc` after delivering the current Euclidean implementation.
Derive a CC-compatible objective/prior scale before implementation; do not reuse
Euclidean whitened-power scaling unchanged. Require independent penalty/derivative
checks, disabled-prior identity, and matched continuous local refinement tests.
This is future scope; the current Euclidean-only contract remains unchanged.

Keep the current width implementation as experimental support; see
[Experimental shift-width policies](#experimental-shift-width-policies) for its
meaning and assumptions. Follow-up work will assess that population model and floor;
defer a production recommendation until multi-dataset validation. Other users may
run datasets unavailable to the developers and report matched prior-off/on results,
width settings/history, sampling/update fractions, starting-pose quality, FSC/cFAR,
convergence and runtime, with source/executable identity and seeds. These reports
provide cross-dataset evidence; availability in the codebase does not establish a
production default or scientific benefit.

#### Corrected background retained for future work

These ideas are deferred, not prerequisites or approved implementation packages:

- **Laplace diagnostics:** Gauss-Newton approximates curvature; it is not the exact
  Hessian. Posterior covariance needs data plus prior curvature in consistent units,
  without LM damping or the calculator's heuristic curvature scaling. Stationarity,
  positive definiteness and active bounds matter. A joint covariance's shift marginal
  is the shift block of the full inverse, not the inverse of its shift block.
  Under these assumptions, `Sigma ~= inverse(T_cart * H_total_raw)`; the calculator's
  unscaled data Gauss-Newton block is `2*H_accum/P`. Local refinement currently has
  no consumer for this covariance, so do not add a getter or persistent output for it.
- **Width feedback:** A per-particle covariance inferred from the same image is not
  automatically an independent prior for the next pass. Keep the population policy
  above. The separate `polar_shift_laplace_width.md` proposal uses widths for search
  support, not as a likelihood term; it is not assumed implemented here.
- **Top-k mixtures:** Refining every candidate and committing the heaviest would
  change the existing random-draw contract. A retained top-k list is not the full
  posterior. Laplace mode weights require comparable objectives, consistent priors,
  valid local covariances and duplicate-mode handling; fitting error alone is not
  mode probability. Cost grows with the number of LM starts and needs measurement.
  The proposed local evidence factor `exp(-E_k)*sqrt(det(Sigma_k))` is meaningful only
  with comparable negative-log posterior energies and the same integration measure;
  separate seed-centered penalties are not automatically one common posterior.
- **Latest separate workflow discussion:** Keep polar extras, optimize only the
  selected pose, and initially retain original polar weights. This supersedes the
  background's all-candidate LM proposal, but is not implemented by this local plan.
- **Future sampling/reweighting:** For inverse-temperature density `exp(-beta*E)`,
  local covariance scales as `Sigma/beta`, not `beta*Sigma`. Normalize future mode
  weights in double-precision log space; this does not add sampling to local LM.

#### Validation gates and deferred decisions

These gates distinguish implementation readiness from later scientific acceptance.

| Item | How to close it | Needed for |
| --- | --- | --- |
| Production width policy | Select from comparative evidence; all-active mean-square is an experimental candidate, not a production decision. Define the optional completed-sample median conversion/population capture before testing it. | Production recommendation; not initial implementation. |
| Numerical tolerances | Developer derives and records tolerances from independent formulas, finite-difference accuracy and pre-change fixtures. Shift-width comparison settings are agreed above; rotation values can wait. | Shift numerical tests. |
| Scientific fixtures and margins | Select calibration/held-out fixtures and budget; derive candidate margins from exploration, then freeze them before independent validation. No additional Hans approval gate. | Scientific validation. |
| Rotation-width interface | Resolved: maintainer approved explicit `cont_prior_rot_sigma`, in per-axis degrees; adaptive rotation remains deferred. | Implemented and tested. |

### 3. Delivery and verification

Launch decision, 10 November 2026: the maintainer chose to finish the guided manual
review and create a reviewed checkpoint before the beta-gal HPC campaign. No dirty
scientific run is authorized. No checkpoint commit or job submission has occurred;
committing still requires the maintainer's explicit instruction after review.

Post-merge Linux evidence supplied 10 November 2026: compilation succeeded;
`local_startup` passed 302 assertions, including both polar/prior rejection cases;
fast gate passed 15/15 and CLI-chain validation (595 keys, 557 bindings, 1,803 UI
inputs). Advisory time was 90.14 s, dominated by unchanged numerics; 86 source
files were reported untracked in the synced Linux checkout. These are build/test
results, not scientific or Git-packaging acceptance.

Beta-gal HPC comparison prepared in the main checkout's local
`.codex/submit_hpc_betagal_priors.sh`; experiment record:
`.codex/diagnostics/cartesian-priors/experiment.md`. Reuses the existing class
checkpoint, regularized particle preparation and one common auto baseline before
five isolated shift-prior arms. No job submitted. Syntax and mock checks passed;
source-review/checkpoint or an explicit noncanonical dirty-run decision precedes
launch. Local-history baseline/final snapshots: `fdaa60c` / `1d46325` (verified).

Retained multi-iteration benchmark verified, 10 November 2026 (pre-upstream-merge
build): `cont_refine3D_bench_mip` completed all six cases in 425.5 s. The configured
budget was six iterations, but every case stopped normally after iteration two
with `CONVERGED: yes` and 98% stable poses. With 100 particles and SNR 0.05,
median rotation/shift errors were 0.2007019 degrees / 0.06283806 A (off),
0.2006799 / 0.06286090 A (adaptive shift), and 0.2007998 / 0.06285888 A (pose).
All cases reported 100% improved poses. Adaptive variance remained at 1.0656
native pixels squared in both iterations: this case does not exercise shrinking
adaptive widths or select a production width policy. No runtime errors were found.
These results do not validate the subsequent upstream integration.

Upstream integration, 10 November 2026: feature worktree fast-forwarded to
`9b2a33167`, matching the main checkout; prior changes reapplied without conflicts.
Upstream adds `refine3D_local refine=prob_neigh` (geom cone). Cartesian priors remain
restricted to `refine=cont`; the local commander now rejects both prior modes on
the polar route before startup. Two rejection cases extend the existing
`local_startup` suite and verify diagnostics and input preservation.
Header, registry, tester-linkage, dependency, import and whitespace checks pass.
Previous Linux results predate this integration; Linux rebuild and regression
checks remain pending. All feature-only files were verified byte-for-byte against
`feature-cartesian-shift-prior-premerge-20261110-092044` before updating this note.

The delivery adds configurable priors. A scientific improvement is not required
for implementation acceptance, and no production width change is proposed.

| Milestone | Evidence |
| --- | --- |
| Design, 6�9 November 2026 | Shift first; joint 5-variable LM; separate Cartesian width; data-only identifiability and reporting; Matrix-Fisher rotation; double-precision prior arithmetic. Rationale is retained in the specification above. |
| Controls approved, 9 November | `cont_prior=none\|shift\|pose`; optional fixed native-pixel shift width; explicit per-axis rotation width in degrees. No inferred rotation default. |
| Independent prior-off baseline | Main-checkout log: 80 checks passed. Six printed records match the feature run exactly: status/counts, rotation, shift, objectives, motion and shell residuals. |
| Linux unit validation | Latest Cartesian suite: 1,310 assertions, eight groups, 7.9 s. Includes 117 optimizer and 86 strategy checks and expected worker/master rejection. |
| Linux fast gate | All 15 tests and CLI-chain checks passed. Total 88.38 s; unchanged numerics suite dominates the advisory overrun. Linux sync reports 85 compiled source files untracked there; this is a packaging caveat. |
| Retained benchmark | Six cases completed in 233.9 s: 100 particles, SNR 0.05, eight threads, one outer iteration, rotation width 5 degrees. Demonstrates execution, not a scientific advantage. |
| Final cleanup | Merged shift/pose transaction setup; corrected displaced width comments; retained exported numerical helpers because private visibility caused observed gfortran/submodule link failures. Focused optimizer rerun passed all 117 assertions (0.2 s), verified from the supplied log. |

Evidence was read from the supplied `unit_cart_align3D.log`, `run_fast_gate.log`,
`pose_optimizer_base.log`, `pose_optimizer_feat.log` and retained benchmark
`priors.log` / `metrics.tsv`. The benchmark median rotation/shift errors were
0.1996218 degrees / 0.04840194 A (off), 0.1995515 / 0.04858338 A (adaptive shift),
and 0.1992852 / 0.04858995 A (pose). These small differences are descriptive.
Adaptive variance was 1.0656 native pixels squared at effective LP 3.019 A;
fixed 1x variance was 1.0519 at requested LP 3 A. They are not identical-width cases.

No Windows build or SIMPLE execution was attempted. Static headers, registry,
tester-linkage, dependency, import and whitespace checks passed before final
cleanup and are rerun for delivery. No commit or push was made.

Before commit: complete the guided manual review. The focused optimizer cleanup
rerun passed. Production shift-width/policy selection remains pending the agreed
test comparison; the smoke benchmark did not choose a winner. Multi-iteration
calibration and independent validation provide that evidence, without requiring
scientific improvement as a prerequisite for prior support. External-scheduler
validation remains separate. Preserve default-off behavior unless separately
supported by evidence and approved.

HPC scheduling (10 November): `submit_hpc_betagal_priors.sh` prepares the common
checkpoint and prior-off control once, then releases a dependent four-task Slurm
array with at most three concurrent arms. Each task uses 8 parts x 12 threads
(96 CPUs) on one node with `SIMPLE_QSYS=local`; peak 288 CPUs.

#### Beta-gal HPC evidence, 10 November 2026

Reviewed downloaded campaign `betagal_priors_20261110_104709_qnC0DR`. Shared
preparation and off control were reused after fixing the launcher stop-marker
check; resumed array tasks all returned zero. All five arms completed six local
iterations and final reconstruction with nonempty final products and no fatal
markers. Resumed task logs show successful input/executable hash checks. The
original preparation exit 2 remains the recorded launcher failure, not a failed
Fortran control run (that command returned zero).

| Arm | Shift width, native pixels | Final FSC 0.143, A | Final cFAR | Local run, seconds |
| --- | --- | --- | --- | --- |
| Off | Disabled | 3.33 | 0.7054 | 542.0 |
| Adaptive | All-active estimator | 3.33 | 0.7050 | 538.9 |
| Fixed 1x | 1 | 3.33 | 0.7050 | 542.1 |
| Fixed 2x | 2 | 3.33 | 0.7055 | 545.3 |
| Fixed 4x | 4 | 3.33 | 0.7051 | 540.5 |

All final FSC 0.5 values are 3.87 A. Every iteration updated all 3,208 active
particles; this does not test stale increments under fractional updates. Stability
was reached from iteration 3 but the workflow enforced a six-iteration minimum.
The saved fixed-width floor is 1 native pixel. Adaptive width history is absent
from the downloaded master/array logs; do not infer that it stayed at the floor.
All arms used the same base-half fallback and sigma bootstrap at startup. The
source was explicitly allowed dirty; Git stamp ee5b68be is not exact binary/source
identity. This is exploratory execution evidence, not independent scientific
acceptance, a width-policy selection, or a speed comparison across cluster loads.

#### Retained adaptive-width diagnostic, 10 November 2026

Each Cartesian prior worker now writes its iteration and actual native-pixel
variance to a small part readout; empty updates write zero. After matching, the
master validates all readouts and prints iteration, nonempty/total part counts
and minimum/maximum sigma in native pixels. This observes worker values without
recomputing the estimator or changing optimization. Missing, stale, unreadable
or nonfinite readouts fail explicitly. Existing local_staging coverage now checks
variance-to-width conversion, empty-part handling and the disabled path.

Launcher `--adaptive-only SAVED_RUN_DIRECTORY` creates a fresh diagnostic output,
reuses the saved preparation through an input link with hash checks, and records new
source/executable hashes for the rebuilt binary. It submits only one 96-CPU local
SIMPLE job; no preparation, control or fixed-width arms run again. Bash syntax,
mock adaptive-only scheduling/checkpoint preservation and normal dependency/array
regression checks passed locally. Fortran compilation, local_staging, fast gate
and the adaptive-only HPC run remain pending on Linux.

#### Upstream integration and review handoff, 10 November 2026

Both local master and feature/cartesian-shift-prior were fast-forwarded from
9b2a33167 to fetched origin/master 749aeccfd (five upstream commits; neither
branch had local commits to replay). Feature changes restored cleanly. Main's
pre-existing generated-index conflicts were resolved by preserving its original
local files; its unrelated work and original three staged plans remain intact.
The independent prior-off fixture in main is unchanged and included once.

The production implementation and essential regression coverage (18 source files)
are now uncommitted in local master for continued manual review. The exploratory
cont_refine3D_bench suite=priors implementation/UI and local HPC helpers remain in
the feature worktree. This single living plan moved to local master with the code;
it is not duplicated. Existing evidence remains here; the optional benchmark is
not part of this delivery's command surface.

Recovery snapshot: sibling cartesian-prior-preintegration-20261110-125158 contains
both original dirty trees, untracked files, staged/worktree patches and hashes.
Static description, test registry, FLEX DAG and whitespace checks passed after
integration. These do not replace Linux compilation or runtime validation of the
merged source, including the new width reporting. HPC compilation is in progress
on the earlier source snapshot. No main-repository commit, push or HPC submission.

#### Delivery scope correction, 10 November 2026

Worker width readouts, report_cartesian_prior_widths, its master call and
local_staging reporting test are diagnostic-only and remain in the feature
worktree. They were removed from local master at the maintainer's request.
The prior calculations and numerical, lifecycle, rollback, propagation and
failure-handling regressions remain in the delivery. The adaptive-only width
experiment must use the feature build; the production build does not emit or
require these readouts. Earlier reporting implementation notes describe that
feature-only diagnostic, not the current master delivery.

#### HPC focused checks reported, 10 November 2026

Reviewed attached console output: local_staging passed 102 assertions (0.76 s);
build_gcc15.2.0 fast gate passed all 15 tests (116.11 s). Command-line chain
passed (595 keys, 557 bindings, 1805 UI inputs). Budget overrun is advisory;
unit_numerics took 116.07 s. The synced checkout warned about 88 untracked compiled
sources. Git stamp ee5b68be does not identify the exact build or establish that
this tested the newly integrated master delivery. These are unit/fast-gate checks;
the adaptive-only real-data width-history result remains pending.

#### Merged-master Linux validation, 10 November 2026

Maintainer confirmed sync and successful Linux compilation of local master after
the transfer and diagnostic removal. Reviewed pasted output: local_startup passed
302 assertions in 6.10 s; fast gate passed all 15 tests in 103.73 s, including
unit_cart_align3D (9.48 s). Command-line chain passed: 595 keys, 557 bindings,
1803 UI inputs. Expected child exit-1 rejection cases did not fail the suite.
The 60-second gate budget warning is advisory (unit_numerics: 103.73 s).
The synced Linux checkout still reports 88 untracked compiled sources; its Git
stamp c97106c4 is not an exact identity for the synced source. Results are tied
to the maintainer's stated master sync and supplied output, not that stamp.
Manual review and feature-only adaptive width-history experiment remain open.
No main-repository commit or push authorized.

Prior-off baseline delivery correction: test_prior_off_baseline and its printed
cross-build records remain only in the feature worktree. They were removed from
master after manual review. The already captured independent baseline evidence
is retained above. test_prior_off_reuse remains an automated regression in master.

#### Focused commit handoff, 10 November 2026

Maintainer completed the guided source review and authorized a local commit,
without push. Local master fast-forwarded to origin/master 4965b3fe8; upstream
even-grid validation and its tests merged without conflicts with the prior work.
The commit contains 17 prior-related Fortran files and this single plan. Local
HPC/width diagnostics, baseline-printing helper and exploratory benchmark are
excluded. The retained benchmark command above is explicitly feature-only.

Static description, test registry, FLEX DAG and whitespace checks passed after
integration. The latest propagation and tight-prior assertions require the
maintainer's post-commit Linux build/tests. Earlier passing runs remain historical
evidence, not a claim that this exact committed snapshot has already passed.
