# One Levenberg-Marquardt for the library

Status, 2026-11-17: planned; revised the same day to one engine for every LM. Base: `e94681175`.
Depends on: nothing for Phases 0-3; the PorF A/B of Phase 4 needs the batch-capacity floor in
`simple_matcher_memory.f90` (working tree, validated: no). Needed by: `polar_continuous_view.md`
(Phases 1-2 here).

## 1. Why

The polar searches refine shifts with L-BFGS-B. In PorF `prob_neigh`, iteration 8, refinement took
950 thread-s against 398 for scoring the bank: 10 shift searches per particle at 13-15 gradient
evaluations each (8.6 in iteration 1). The Cartesian route's LM needs 3.3 calculator passes per
five-parameter solve by iteration 8 (9.0 in iteration 1). The factory offers `lbfgsb`, `simplex` and
`de`. The library has three LM loops, none in the factory: `refine_joint` and `refine_shift` in
`cartft_pose_opt`, and the Guinier fit of `simple_dose_series`. Rulings (Hans, 2026-11-17):

- both solvers everywhere a client can use them;
- one LM in the library;
- each client's measurement decides its default.

**Checked at the base.**

- **Factory clients with derivatives:** `pftc_shsrch_grad` (shift mode, 2 variables; joint mode, 3),
  `ftexp_shsrch` (motion correction), `ctf_estimate_cost4Dcont`, `quality_logistic`.
  `ctf_estimate_cost2D` (DE) and `volpft_symsrch` (simplex) use no derivatives.
- **Polar Euclidean score:** exp(−L), and the shift prior multiplies it by exp(−w|d|²), so
  L + w|d|² has the minimiser of today's cost. The cc shift costs divide by norms that do not
  depend on the shift.
- **`refine_joint` is generic except for four pose-specific steps:**
  - the SO(3) right increment;
  - the rotation and shift step caps;
  - the cumulative guard against the transaction seed;
  - the prior, added after the data-only identifiability test.
- **`refine_shift`** (a 2×2 closed-form solve) serves only the internal `cont_route=shift_then_joint`
  (production runs `joint`) and the tester's shift-only cases.
- **The Guinier fit** is classic Marquardt (λ ×10 on a failed try, ×0.1 after a success), with 2
  unbounded parameters.
- **`simple_linalg`** owns `cholesky`. Its loop order is that of `solve_pose_cholesky`.

## 2. Contract

- **One LM algorithm.** After Phase 3 no other LM loop exists in the library.
- **The engine owns the algorithm only.** Each client supplies its terms, its step limits, its trial
  point and its constraints.
- **The Cartesian route is byte-identical** at one thread after Phase 1. The exception is the
  shift-only path (decision 10).
- **Factory protocol unchanged.** `specify('lm', ...)`, `set_fdfcurvfun_8`, factory `new`, `x` and
  limits, `minimize`. The result is in `x`/`x_8`, the counts in `nevals`, `ngevals`, `niter`.
- **No default changes before Phase 4's read.** The dose-series fit is the one result-changing
  migration, judged by its own suites.
- One problem or spec per thread, as now.

## 3. Decisions (2026-11-17)

1. **Two layers, following `nlcg_solver`.**
   - **The engine,** `src/utils/math/simple_lm_solver.f90`. `lm_solve(problem, options, outcome)`
     works on an abstract `lm_problem` whose deferred procedures belong to the client:
     - `eval`: f, gradient and curvature at the current trial point;
     - `limit_step`: may shorten a step, and reports whether it did;
     - `propose`: forms the trial point and reports whether it is feasible;
     - `accept`: makes the trial current, receives the ratio, and reports whether its objective
       changed.

     `penalty` (prior terms, added after the identifiability test) defaults to nothing. The options
     are `maxits`, a coordinate `scale(D)` and `active(D)`. Optional seed terms skip the first
     evaluation. The outcome holds the status and the counters. The engine lives in `src/utils` so
     that `simple_dose_series` can import it too.
   - **The factory optimizer,** `opt/simple_opt_lm.f90`. `opt_lm` extends `optimizer`, is selected by
     `'lm'` and works in double precision. It runs a private `lm_problem` over `opt_spec`: `eval` calls
     `fdfcurvfun_8`, `limit_step` clips to `limits`, `propose` adds the step, and `accept` copies it
     and runs `opt_callback`.
2. **Callback.** `opt_spec` gains `fdfcurvfun_8(fun_self, vec, f, grad, curv, D)` and
   `set_fdfcurvfun_8`: f, the gradient and a symmetric positive semi-definite curvature. `'lm'` stops
   at `new` without it. A client that sets both callbacks runs under either solver.
3. **Step control: `refine_joint`'s, extracted verbatim.**
   - **Damped system.** C + μ·diag(d), with d_i = max(C_ii, √ε·max|C|, ε²), solved by `cholesky`
     with pivot floor √ε·max|·|.
   - **Damping updates.** μ₀ = 1e-3. A rejection multiplies μ by m, where m starts at 4 and doubles;
     an acceptance resets m. μ halves when the ratio exceeds 0.75.
   - **Acceptance.** A step is accepted when the actual reduction is positive and at least 0.25 of
     the predicted.
   - **Stops.** 8 consecutive rejections; predicted reduction below 1e-6·|f|; accepted relative
     reduction below 1e-6; |g| below 1e-8; `maxits` iterations.
   - **Outcomes.** Its six outcomes (accepted, no improvement, no reliable update, bound rejected,
     invalid numerics, iteration limit), to which `CARTFT_*` map one-to-one.
4. **Bounds (factory).** `limits` are required. `limit_step` clips the step to the box, and the
   predicted reduction is the clipped step's, so a minimum on a bound is reached as `lbfgsb` reaches
   it.
5. **Hook (factory).** `opt_callback` runs in `accept`. When it changes the objective (the
   alternating search's discrete angle), the engine re-evaluates and counts the evaluation.
6. **Curvature by objective.**
   - **Euclidean.** The whitened residual's Gauss-Newton matrix JᵀWJ. The shift prior is a
     `penalty` with curvature 2w·I.
   - **cc with shift-invariant norms** (polar shift search at a fixed rotation, `ftexp_shsrch`). The
     Gauss-Newton matrix of the equivalent least-squares problem. It is constant over the shift, so it
     is computed once per solve.
   - **cc whose norms move** (the polar joint search under cc, the CTF fit). The normalised model's
     Gauss-Newton matrix, following `cartft_calc`'s `finalize_cc_normal_terms`.
   - **Log-likelihood** (`quality_logistic`). The Fisher matrix plus the ridge, which is its exact
     Hessian.
7. **Per-client switch.** Each factory client's solver becomes a named private constant beside its
   `specify`, and its `new` takes an optional solver name that only testers pass. There is no CLI key:
   real-data A/Bs use two builds from one snapshot.
8. **Per-client verdict, pre-registered.** LM wins on a client's fixture if all of these hold:
   - median thread-seconds per solve at most 0.8× L-BFGS-B's (guessed);
   - accuracy within the client's existing floor and within the base's spread of L-BFGS-B's
     (measured in Phase 0);
   - no failure that L-BFGS-B does not have.

   Otherwise the client keeps L-BFGS-B, and the negative is recorded with its mechanism.
9. **The Cartesian route on the engine.** `refine_joint` becomes `lm_solve` on a
   `cartft_pose_problem` in `cftc/`:
   - `eval` is `pose_normal_terms`, keeping the trial's residual;
   - `penalty` is `pose_prior_terms`;
   - `limit_step` applies the rotation and shift norm caps;
   - `propose` applies the right increment and the cumulative guard;
   - `accept` adopts the pose, its terms and its residual, and bins the ratio.

   The scale is (rs, rs, rs, 1, 1) and `active` follows `trs`. The transaction in `refine_pose` (seed
   evaluation, final acceptance on the penalised loss) is unchanged. The LM helpers of
   `cartft_pose_opt` and `_system` move into the engine: `build_pose_lm_system`,
   `solve_pose_cholesky`, `increase_lm_damping`, and the scaling and masking.
10. **`shift_then_joint` retired; a shift-only path kept (ruled).**
    - **Retired.** The internal key `cont_route`, `l_cont_shift_first`, `new`'s `shift_first` and
      the shift-first branch of `refine_pose` go. The `simple-dead-code` compiler probe lists every
      reference to remove.
    - **Kept.** A shift-only Cartesian solve stays for Cyril's motion-polishing code. `refine_shift`
      becomes public with its current signature (fixed rotation, shift in/out), implemented as
      `lm_solve` on a 2-variable `cartft_shift_problem`:
      - `eval` is `shift_normal_terms`;
      - `limit_step` applies the shift step cap;
      - `propose` applies the cumulative shift bound.

      Its numbers change (the engine's Cholesky and tests replace the closed-form 2×2 solve). With
      no production caller, it is judged by the tester's existing shift-only recovery cases at their
      tolerances.
11. **Dose series.** The Guinier fit becomes an `lm_problem` (2 parameters, unbounded). This changes
    results, judged by `unit_numerics` suites `dose series` and `dose series manifest` at their
    existing floors.

Defaults taken without a question: the names `lm_problem`, `lm_solve`, `opt_lm`, `fdfcurvfun_8`.

**Not in this plan:** the DE and simplex clients; the NLCG and PCG engines; the hybrid's view
refinement (`polar_continuous_view.md`).

## 4. Acceptance and stops

The LM tests go in `simple_opt_tester.f90` (`unit_numerics suite=optimisers`), through the factory:

- **L1.** Rosenbrock in residual form, r = (10(y − x²), 1 − x), from the tester's start, reaches
  (1, 1) within 1e-6 (guessed). Evaluations against `lbfgsb` are reported, not gated.
- **L2.** A quadratic whose unconstrained minimum lies outside the box ends on the bound, equal to
  `lbfgsb`'s within 1e-6 (guessed).
- **L3, negative control.** With the gradient's sign flipped, LM accepts no step and returns the start
  after 8 rejections.
- **L4.** `opt_callback` runs once per accepted step, and `nevals` counts the re-evaluation.

**Curvature oracle, per client.** At zero residual the curvature equals the central-difference
Hessian of the analytic gradient, within the client's gradient tolerance, and it is positive
semi-definite. The logistic Fisher matrix must equal the Hessian at any point.

**Stops.**

- The Cartesian identity fails: diagnose, never accept a near-identity, and ask if the cause is the
  `cholesky` swap.
- A client's oracle fails: exclude that client and report it.
- An LM test passes only with a constant other than `refine_joint`'s.
- A file outside the table.

**Gates no ruling waives.** The fast gate; the optimisers suite with L3; the Cartesian identity
(Phase 1); the polar identity (Phase 3); `polar_matching_bench suite=fills`.

## 5. Phases

| Phase | Kind | Content | Exit |
| --- | --- | --- | --- |
| 0 | measurement | Release + Debug builds; fast gate; `unit_numerics`, `unit_image`, `unit_pftc_align2D3D`, `unit_cart_align3D`, `lib_cart_align3D`, `fills`. Identity cases at `nthr=1`, run twice: `refine3D_local refine=cont` and `refine=prob_neigh` on the `cont_refine3D_1jxy` project. Per-client baselines (evaluations, thread-seconds, accuracy). The dead-code probe's reference list for `cont_route`. Hans names the motion-correction and CTF sets | Recorded in Progress |
| 1 | identical (the shift-only path aside) | Engine extracted from `refine_joint`; the Cartesian route on it; `shift_then_joint` retired; `refine_shift` public on the engine | Cartesian identity byte-identical; cart suites print identical numbers, the shift-only cases pass at their tolerances; `rg` finds none of the moved helpers in `cftc/` and no `cont_route` in `src/` |
| 2 | tests-only | `opt_lm`, the spec callback, the factory entry | L1-L4 pass; `'lm'` has no production caller |
| 3 | tests-only (factory clients, production identical); result-changing (dose series) | Curvatures and the solver name for the five factory clients, each tester running both solvers; the Guinier fit on the engine | Oracles pass; one comparison table per client; polar identity byte-identical; dose-series suites at their floors |
| 4 | result-changing | Decision 8 per client. For each winner: the default switched in build B; its workflow gate (`simulated_workflow` profile for the shift searches); a real-data A/B (PorF `local4_best` command, two runs per build, or the named sets); Hans reads, and a default switches only on his ruling. Report `lm_optimizer_report.md` | Verdicts, A/B numbers and rulings in Progress |

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/utils/math/simple_lm_solver.f90` (new) | The engine | 1 |
| `cftc/simple_cartft_pose_opt.f90`, `simple_cartft_pose_opt_system.f90`, `simple_cartft_pose_opt_tester.f90` | Cartesian pose and shift problems; helpers moved out; public `refine_shift` | 1 |
| `params/simple_parameters{,_parse,_phases,_tester}.f90`, `strategies/search/simple_strategy3D_cont{,_1jyx_tester}.f90`, `commanders/test/simple_commanders_test_highlevel_cont.f90`, `ui/simple_ui_visibility_tester.f90` (if needed) | Retire `cont_route` | 1 |
| `opt/simple_opt_lm.f90` (new), `simple_opt_spec.f90`, `simple_opt_factory.f90`, `simple_opt_tester.f90` | Factory LM, callback, L1-L4 | 2 |
| `pftc/simple_pftc_shsrch_grad.f90`, `simple_polarft_calc.f90`, `simple_polarft_corr.f90`, `simple_pftc_inplane_tester.f90` | Shift and joint curvature, solver name | 3-4 |
| `image/simple_ftexp_shsrch.f90`, `_tester`, `simple_ft_expanded.f90` (if needed); `ctf/simple_ctf_estimate_cost.f90`, `simple_ctf_tester.f90`; `quality/simple_quality_logistic.f90` | Curvatures, solver name | 3-4 |
| `src/utils/filter/simple_dose_series.f90` | Guinier fit on the engine | 3 |
| `doc/policies/3D/refine3D_policy.md`; skills `simple-main-opt`, `simple-main-cftc`, `simple-refine3d` (approved with this plan) | Living descriptions | 1-2 |
| this note | Progress | 0-4 |

Paths under `src/main/` unless rooted.

## Progress

- 2026-11-17, session: plan written; premises checked at `e94681175`; nothing built. Revised the same
  day: one engine below every LM client, the Cartesian route and the dose-series fit included.
