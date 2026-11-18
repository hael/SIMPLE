# One Levenberg-Marquardt for the library

Status, 2026-11-18: Phases 0-3 done (run `lm_optimizer`; interim report `lm_optimizer_report.md`); Phase 4 open.
Planned 2026-11-17; revised the same day to one engine for every LM, then after the design
review (section 6). Base: `aae467bc7` (numerical source unchanged from `e94681175`).
Depends on: nothing for Phases 0-3; the PorF A/B of Phase 4 needs the batch-capacity floor
(`aae467bc7`, validated: no). Needed by: `polar_continuous_view.md` (Phases 1-2 here).

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
- **`simple_linalg`** owns `cholesky`, `chol_forward` and `chol_backward` (users: flex and XD-GMM).
  `cholesky` subtracts each product in turn, while `solve_pose_cholesky` subtracts a `dot_product`.
  The two round differently (review 3), so swapping one for the other does not preserve identity.

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

     `penalty` (prior terms, added after the observability test) defaults to nothing. Two hooks
     default to "always true" (revised after review 1 and 2):
     - `observable` receives the scaled data system and may veto the solve;
     - `free` marks the coordinates the solve may move this iteration.

     The engine solves on `active .and. free`. The options are `maxits`, a coordinate `scale(D)` and
     `active(D)`. Optional seed terms skip the first evaluation. The outcome holds the status and the
     counters. The engine lives in `src/utils` so that `simple_dose_series` can import it too.
   - **The factory optimizer,** `opt/simple_opt_lm.f90`. `opt_lm` extends `optimizer`, is selected by
     `'lm'` and works in double precision. It runs a private `lm_problem` over `opt_spec`:
     - `eval` calls `fdfcurvfun_8`;
     - `free` freezes the coordinates that are active on a bound (decision 4);
     - `limit_step` clips to `limits`;
     - `propose` adds the step;
     - `accept` copies it and runs `opt_callback`.
2. **Callback.** `opt_spec` gains `fdfcurvfun_8(fun_self, vec, f, grad, curv, D)` and
   `set_fdfcurvfun_8`: f, the gradient and a symmetric positive semi-definite curvature. `'lm'` stops
   at `new` without it. A client that sets both callbacks runs under either solver.
3. **Step control: `refine_joint`'s, extracted verbatim** (revised after review 2 and 3).
   - **Damped system.** C + μ·diag(d), with d_i = max(C_ii, √ε·max|C|, ε²). It is solved by a new
     `simple_linalg` routine, `spd_solve_rel`: `solve_pose_cholesky` generalised to n, with its
     reduction order, its relative pivot floor √ε·max|·| and its finite checks.
   - **Damping updates.** μ₀ = 1e-3. A rejection multiplies μ by m, where m starts at 4 and doubles;
     an acceptance resets m. μ halves when the ratio exceeds 0.75.
   - **Acceptance.** A step is accepted when the actual reduction is positive and at least 0.25 of
     the predicted.
   - **Unreliable solves.** A damped solve that fails the pivot test counts as a rejection (μ grows)
     instead of ending the solve. Today `refine_joint` stops; PorF showed 0 no-reliable-update
     outcomes in about 300,000 solves.
   - **Stops.** 8 consecutive rejections; predicted reduction below 1e-6·|f|; accepted relative
     reduction below 1e-6; masked scaled |g| below 1e-8; `maxits` iterations. (Amended by R3: the three
     tolerances are `lm_options`, these values by default; the factory passes its spec's.)
   - **Outcomes.** Its six outcomes (accepted, no improvement, no reliable update, bound rejected,
     invalid numerics, iteration limit), to which `CARTFT_*` map one-to-one.
4. **Bounds (factory): a solve on the free coordinates** (revised after review 1). Componentwise
   clipping of a coupled step can point uphill: the review's counterexample stalls at the start.
   - **Free coordinates.** Each iteration, the adapter's `free` freezes a coordinate that sits on a
     bound with its gradient pointing outward.
   - **Clipping.** `limit_step` clips the free step to the box, and the predicted reduction is the
     clipped step's.
   - **Convergence.** The stationarity test on the masked gradient is then the projected-gradient
     test. As μ grows, the clipped step tends to scaled steepest descent on the free coordinates,
     which descends, so a minimum on a bound is reached.
5. **Hook (factory).** `opt_callback` runs in `accept`. When it changes the objective (the
   alternating search's discrete angle), the engine re-evaluates and counts the evaluation.
6. **Curvature by objective** (revised after review 2 and 4).
   - **Euclidean.** The whitened residual's Gauss-Newton matrix JᵀWJ.
   - **Polar shift prior.** In the shift mode (a fixed rotation) it is a `penalty` with curvature
     2w·I. In the joint mode the centre c(a) = R(a)u0 turns with the angle, so the curvature is
     2w·[I, −c′(a)]ᵀ[I, −c′(a)], in rotation-index units (`prior_centre_angle`'s dc0). That is
     also the exact Hessian at zero prior residual.
   - **cc with shift-invariant norms** (polar shift search at a fixed rotation, `ftexp_shsrch`). The
     Gauss-Newton matrix of the equivalent least-squares problem. It is constant over the shift, so it
     is computed once per solve.
   - **cc whose norms move** (the polar joint search under cc, the CTF fit). The normalised model's
     Gauss-Newton matrix, following `cartft_calc`'s `finalize_cc_normal_terms`.
   - **Log-likelihood** (`quality_logistic`). The Fisher matrix plus the ridge, which is its exact
     Hessian.
   - **Null directions** (replaced by R4). A positive semi-definite curvature may be singular: in the
     CTF fit at dfx = dfy the astigmatism angle has a zero Jacobian. After the observability test, a
     coordinate whose scaled gradient and scaled curvature diagonal are both exactly zero gets the
     inactive mask's unit diagonal in the damped solve: its step is zero and the others move. (The
     damping floor alone does not do it: such a coordinate never passes the relative pivot floor, Phase 1
     finding.) The generic engine vetoes nothing.
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
   - `observable` is the data-only Cholesky test of `prepare_pose_data_system` on the scaled, masked
     data system, so a prior can never rescue an unobservable pose (revised after review 2);
   - `limit_step` applies the rotation and shift norm caps;
   - `propose` applies the right increment and the cumulative guard;
   - `accept` adopts the pose, its terms and its residual, and bins the ratio.

   The scale is (rs, rs, rs, 1, 1) and `active` follows `trs`. The transaction in `refine_pose` (seed
   evaluation, final acceptance on the penalised loss) is unchanged. `build_pose_lm_system`,
   `increase_lm_damping`, and the scaling and masking move into the engine; `solve_pose_cholesky`
   becomes `spd_solve_rel` in `simple_linalg`. `test_prior_identifiability` and
   `test_invalid_and_unobservable_inputs` stay as they are.
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

Defaults taken without a question: the names `lm_problem`, `lm_solve`, `opt_lm`, `fdfcurvfun_8`,
`spd_solve_rel`, `observable`, `free`.

Rulings (Hans, 2026-11-17; first applied in Phase 0):

- **R1. Fixtures only in Phase 0.** The Phase 0 baselines for motion correction (`ftexp_shsrch`) and
  the CTF fit (`ctf_estimate_cost4Dcont`) use the existing test fixtures only. The real
  motion-correction and CTF sets move to Phase 4, where Hans names them.
- **R2. This run is Phases 0-3.** Phase 4 (decision 8's verdicts, the default switches, the real-data
  A/B and their rulings) is a later run with Hans. The plan stays in `planned/`, and the report of
  Phases 0-3 is an interim report that Phase 4 extends.
- **R3 (2026-11-18; first applied in Phase 2). The engine's stopping tolerances are options.**
  `lm_options` gains the predicted-reduction stop (relative to |f|), the accepted relative-reduction stop
  and the masked scaled gradient stop, defaulting to `refine_joint`'s 1e-6, 1e-6 and 1e-8, which the
  Cartesian problems keep. `opt_lm` takes them from `opt_spec` with L-BFGS-B's meaning: both relative
  stops at `factr`·ε(double), the gradient stop at `pgtol` on the projected (free) gradient. L2 and L3
  set `factr` and `pgtol` explicitly and keep their criteria. Amends decision 3; the section 4 stop on
  `refine_joint`'s constants applies to the Cartesian route only.
- **R4 (2026-11-18; first applied in Phase 2). Null coordinates get a zero step.** Replaces decision 6's
  "Null directions" bullet; new test L6. Both are implemented in Phase 2 (the engine and `opt_lm`), and
  Phase 1's Cartesian identity gate is rerun before Phase 2 is done.

**Not in this plan:** the DE and simplex clients; the NLCG and PCG engines; the hybrid's view
refinement (`polar_continuous_view.md`); merging `spd_solve_rel` with `cholesky`, which would move
the flex and XD-GMM numbers.

## 4. Acceptance and stops

The LM tests go in `simple_opt_tester.f90` (`unit_numerics suite=optimisers`), through the factory:

- **L1.** Rosenbrock in residual form, r = (10(y − x²), 1 − x), from the tester's start, reaches
  (1, 1) within 1e-6 (guessed). Evaluations against `lbfgsb` are reported, not gated.
- **L2, coupled bound (review 1).** f = ((x + 2y + 1)² + (y − 1.9)²)/2 on [0, 10] × [−10, 10],
  starting at (0, 0) on the bound, reaches (0, −0.02) within 1e-6 (guessed), and its projected gradient
  falls below 1e-8. Plain clipping must fail this case: that is the negative control.
- **L3, negative control.** With the gradient's sign flipped, LM accepts no step and returns the start
  after 8 rejections.
- **L4.** `opt_callback` runs once per accepted step, and `nevals` counts the re-evaluation.
- **L5, rank deficiency (review 2).** r = x + y − 1 with a singular curvature reaches f < 1e-12, and
  x − y keeps its start value: the null direction gets no step.
- **L6, null coordinate (R4).** r = x − 1 in two variables reaches f < 1e-12, y keeps its start value,
  and the status is not "no improvement".

**Curvature oracle, per client.** At zero residual the curvature equals the central-difference
Hessian of the analytic gradient, within the client's gradient tolerance, and it is positive
semi-definite. The logistic Fisher matrix must equal the Hessian at any point. The polar joint
oracle runs with a nonzero prior centre and the angle active. The CTF client also runs from an
isotropic start (dfx = dfy): the defocus improves and the angle does not move.

**Stops.**

- The Cartesian identity fails: diagnose it, and never accept a near-identity.
- A client's oracle fails: exclude that client and report it.
- An LM test passes only with a constant other than `refine_joint`'s (since R3: on the Cartesian
  route only; the factory's tests set `factr` and `pgtol`).
- A file outside the table.

**Gates no ruling waives.** The fast gate; the optimisers suite with L3; the Cartesian identity
(Phase 1); the polar identity (Phase 3); `polar_matching_bench suite=fills`.

## 5. Phases

| Phase | Kind | Content | Exit |
| --- | --- | --- | --- |
| 0 | measurement | Release + Debug builds; fast gate; `unit_numerics`, `unit_image`, `unit_pftc_align2D3D`, `unit_cart_align3D`, `lib_cart_align3D`, `fills`. Identity cases at `nthr=1`, run twice: `refine3D_local refine=cont` and `refine=prob_neigh` on the `cont_refine3D_1jxy` project. Per-client baselines (evaluations, thread-seconds, accuracy). Cartesian no-reliable-update counts on the identity case and the cart suites, split between the observability veto and a failed damped solve (any damped-solve failure: stop and ask, decision 3). The dead-code probe's reference list for `cont_route`. The motion-correction and CTF baselines on the existing test fixtures (R1) | Recorded in Progress |
| 1 | identical (the shift-only path aside) | `spd_solve_rel` in `simple_linalg`; the engine extracted from `refine_joint`; the Cartesian route on it; `shift_then_joint` retired; `refine_shift` public on the engine | Cartesian identity byte-identical; the cart suites print identical numbers, and the prior and unobservable tests pass unchanged; the shift-only cases pass at their tolerances; `rg` finds none of the moved helpers in `cftc/` and no `cont_route` in `src/` |
| 2 | tests-only | `opt_lm` with the `free` hook, the spec callback, the factory entry; R3 and R4 in the engine | L1-L6 pass; `'lm'` has no production caller; the Cartesian identity rerun byte-identical (R3, R4) |
| 3 | tests-only (factory clients, production identical); result-changing (dose series) | Curvatures and the solver name for the five factory clients, each tester running both solvers; the Guinier fit on the engine | Oracles pass; one comparison table per client; polar identity byte-identical; dose-series suites at their floors |
| 4 | result-changing | Decision 8 per client. For each winner: the default switched in build B; its workflow gate (`simulated_workflow` profile for the shift searches); a real-data A/B (PorF `local4_best` command, two runs per build, or the named sets); Hans reads, and a default switches only on his ruling. Report `lm_optimizer_report.md` | Verdicts, A/B numbers and rulings in Progress |

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/utils/math/simple_lm_solver.f90` (new) | The engine; R3 tolerances and R4 null coordinates | 1, 2 |
| `src/utils/math/simple_linalg.f90`, `simple_linalg_tester.f90` | `spd_solve_rel`, moved from `cftc` and generalised to n | 1 |
| `cftc/simple_cartft_pose_opt.f90`, `simple_cartft_pose_opt_system.f90`, `simple_cartft_pose_opt_tester.f90` | Cartesian pose and shift problems; helpers moved out; public `refine_shift` | 1 |
| `params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90`, `simple_parameters_tester.f90`, `strategies/search/simple_strategy3D_cont.f90`, `simple_strategy3D_cont_1jyx_tester.f90`, `commanders/test/simple_commanders_test_highlevel_cont.f90`, `ui/simple_ui_visibility_tester.f90` (if needed) | Retire `cont_route` | 1 |
| `opt/simple_opt_lm.f90` (new), `simple_opt_spec.f90`, `simple_opt_factory.f90`, `simple_opt_tester.f90` | Factory LM, `free` hook, callback, L1-L5 | 2 |
| `pftc/simple_pftc_shsrch_grad.f90`, `simple_polarft_calc.f90`, `simple_polarft_corr.f90`, `simple_pftc_inplane_tester.f90` | Shift and joint curvature, solver name | 3, 4 |
| `image/simple_ftexp_shsrch.f90` (and its `_tester`), `simple_ft_expanded.f90` (if needed); `ctf/simple_ctf_estimate_cost.f90`, `simple_ctf_tester.f90`; `quality/simple_quality_logistic.f90` (and its `_tester`) | Curvatures, solver name | 3, 4 |
| `src/utils/filter/simple_dose_series.f90` (and its `_tester`) | Guinier fit on the engine | 3 |
| `doc/policies/3D/refine3D_policy.md`; skills `.github/skills/simple-main-opt/SKILL.md`, `.github/skills/simple-main-cftc/SKILL.md`, `.github/skills/simple-refine3d/SKILL.md` (approved with this plan) | Living descriptions | 1, 2 |
| `doc/implementation_notes/planned/lm_optimizer.md` (this note), `lm_optimizer_report.md` | Progress; the report | 0, 1, 2, 3, 4 |

Paths under `src/main/` unless rooted. A listed `x.f90` covers its `x_tester.f90` (the run's fence
reads the table that way). Phases are listed one by one, not as ranges, for the same reason.

## 6. Design review (2026-11-17)

Reviewed against the current source at `aae467bc7`; the relevant numerical source is unchanged
from `e94681175`. The engine/factory split is appropriate, but the following contracts need
revision before implementation. These findings do not change the decisions above.

1. **P1 — Box clipping is not sufficient for the promised bound convergence (decision 4).**
   A coupled unconstrained step can point uphill after componentwise clipping. Independent
   double-precision Python emulation of the stated controls gives a concrete counterexample:
   `f(x,y) = ((x + 2y + 1)^2 + (y - 1.9)^2)/2`, box `[0,10] x [-10,10]`, start `(0,0)`.
   Its curvature is positive definite, `[[1,2],[2,5]]`, but the solve stops at iteration 6 on
   predicted reduction with no accepted step. The exact constrained minimum is `(0,-0.02)`;
   the feasible y-gradient at the returned start is 0.1. Specify a solve on the free coordinates
   or a projected descent fallback and a projected-gradient convergence check. L2 must include
   a coupled quadratic starting on a bound, not just a diagonal quadratic.

2. **P1 — Data identifiability is a Cartesian policy, not a generic LM prerequisite (decisions 1–3).**
   `prepare_pose_data_system` Cholesky-tests the undamped data matrix and rejects a singular
   active block before adding either a prior or damping. Extracting that rule verbatim conflicts
   with the callback's positive-semidefinite contract. In the supported 3/4-variable CTF fit,
   `dfx=dfy` makes the astigmatism-angle Jacobian identically zero
   (`simple_ctf_estimate_cost::fdf`), although the defocus coordinates can still improve.
   Keep the Cartesian observability veto in its adapter, and define how the generic engine
   regularizes or freezes null directions, including its relative pivot threshold. Add a
   rank-deficient least-squares check and an isotropic-CTF start; preserve the existing Cartesian
   test that forbids a prior from rescuing unobservable pose data.

3. **P2 — The Cholesky swap does not preserve the stated arithmetic order (premise and decision 3).**
   `solve_pose_cholesky` traverses rows and computes `A(i,j) - dot_product(...)`;
   `simple_linalg::cholesky` traverses columns and repeatedly subtracts each product from
   `A(i,j)`. They also differ in finite-value checks. A fixed-seed 5x5 SPD Python comparison
   of those arithmetic forms produced a differing factor entry (5.55e-17 absolute).
   This does not establish a compiled workflow failure, but it disproves the premise used
   to justify byte identity. Preserve the Cartesian reduction order and checks in the shared
   numerical owner, with the necessary file-table change, or separate this substitution from
   the identity-preserving extraction. Do not weaken the identity gate.

4. **P2 — The polar joint prior needs angle and cross curvature (decision 6).**
   `joint_grad_at_angle` uses `P = w |s - c(a)|^2`, where `c(a)=R(a)u0` follows the
   stored seed. For a nonzero seed, its PSD Gauss-Newton curvature is
   `2w [I,-c'(a)]^T [I,-c'(a)]`, including shift-angle and angle-angle entries; `2w I`
   describes only fixed-angle shifts. At zero prior residual this full matrix is also the
   exact Hessian. Specify it in the client contract, using the existing rotation-index units,
   and extend the curvature oracle to a nonzero prior centre with angle active. Otherwise the
   proposed matrix cannot pass the stated zero-residual Hessian oracle for that client.

Validation: source inspection and independent standard-library Python numerical checks only.
No Fortran implementation, build, fast gate, workflow identity, real-data or Linux/BOX run was
performed for this design review. The revised contracts and their compiled regressions remain open.

**Disposition (Hans, 2026-11-17: "update the plan accordingly").** All four findings were accepted.

1. Decision 4 now solves on the free coordinates (the `free` hook); L2 is the review's coupled case.
2. The observability veto moves into the Cartesian problem (`observable`, decisions 1 and 9). The
   generic engine handles null directions by the damping floor, and an unreliable solve raises μ
   (decisions 3 and 6). L5, the isotropic CTF start and the two Cartesian prior/unobservable tests
   check this.
3. The Cartesian factorization moves to `simple_linalg` as `spd_solve_rel`, with its reduction order
   and checks, and the identity gate is unchanged (decision 3; premise corrected). Merging it with
   `cholesky` is out of scope.
4. Decision 6 gives the joint-mode prior curvature 2w·[I, −c′]ᵀ[I, −c′], and the polar joint oracle
   runs with a nonzero centre and the angle active.

## Progress

- 2026-11-17, session: plan written; premises checked at `e94681175`; nothing built. Revised the same
  day: one engine below every LM client, the Cartesian route and the dose-series fit included.
- 2026-11-17, session: design review (section 6) accepted. Decisions 1, 3, 4, 6 and 9, the premise
  on `simple_linalg`, the tests (L2, L5, oracles), Phases 0-2 and the file table were revised; base
  re-pinned to `aae467bc7`. Nothing built.
- **2026-11-17, Phase 0 (measurement), arun `lm_optimizer`.** FWM-S116118: macOS 26.7.2, Apple M2 Ultra,
  24 cores, 192 GiB; Homebrew GCC 16.2.0 (`gfortran-16`, `gcc-16`, `g++-16`, chosen by
  `cmake/AutoDetectGCC.cmake`), cmake 4.2.1; builds at `-j12`. Base `aae467bc7` exported with
  `git archive` to `scratch/keep/base_src`. Load at start 4.7; the process lists taken during the session
  showed no one else's compute jobs, only OS daemons (indexing, antivirus); free disk 5.9 TiB. Paths below are under `/Users/elmlundho/agent_work/arun/lm_optimizer`.
  Files changed: this note only (R1 and R2 recorded in section 3; the Phase 0 row reads R1).
  - **Builds and gate.** Release and Debug base trees `scratch/keep/base_release`, `scratch/keep/base_debug`
    (`scratch/logs/build_base_{release,debug}.log`; no compiler warning). Fast gate 15/15 on both:
    Release 18.7 s, Debug 54.4 s (`scratch/logs/fastgate_base_{release,debug}.log`).
  - **Suites (Release, 1 thread, CTest environment; `scratch/runs/suites_rel1/`).** `unit_numerics` 31/31
    groups, 2060 assertions; `unit_image` 14/14, 976; `unit_pftc_align2D3D` 15/15, 6327; `unit_cart_align3D`
    8/8, 1183; `lib_cart_align3D` (8 threads) 106/106 (joint 54, `shift_then_joint` 52), 24 s;
    `unit_reconstruction` (fast gate; 1378 assertions, `dose series manifest` 87, in the instrumented runs
    below); `polar_matching_bench suite=fills nthr=8`:
    48 rows, no differing entry in any column, 329 s (`scratch/logs/fills_base_release.tsv`).
  - **Identity cases.** The test keeps its seeded project in memory and deletes its fixture on PASS, so
    `scratch/capture_fixture.sh` ran `simple_test_exec test=cont_refine3D_1jxy nthr=1` (the commander honours
    `nthr`) on the Release base, kept every version of stage a's project and stopped the test at its first
    iteration: version 1 is `enter_stage`'s copy of the seeded project (later versions add `eo_byindex` and
    the sigma2 registration). Fixture: `scratch/keep/identity_fixture/seeded_cont_gate.simple` (md5
    `042c7c1b…`) with the simulated stack and maps of `test_cont_refine3D_1jxy_15498/`. Each case,
    `scratch/run_identity.sh <tree> <dir> <refine>`, is stage e's command with `refine=` added:
    `refine3D_local vol1=1JYX.mrc mkdir=no pgrp=c1 mskdiam=120 nsample=500 maxits=4 nthr=1 nparts=1`,
    `SIMPLE_SEED=20260923` (the test seed), a fresh directory (`nparts`, required by the UI, makes it a
    one-part distributed run). Run twice each: `refine=cont` 63 s, `refine=prob_neigh` 290 s.
    `scratch/compare_runs.sh` (`scratch/logs/cmp_base_{cont,prob_neigh}.txt`): in both cases 63 files are
    byte-identical between the two runs, among them every `.mrc`, `fsc*.bin`, `sigma2_state.bin`,
    `simple_stats.txt`, `RESOLUTION_*`, `trailrec_state01.txt`, the Cartesian or reprojection reference
    files, and the text dumps of the `ptcl3D` and `stk` segments (`print_project_field`).
    **Identity outputs:** those byte-identical files. **Declared exclusions** (fixed here, never extended):
    `cont_gate.simple` itself (absolute paths; judged by its `ptcl3D`/`stk` dumps, bytes, and its `out` and
    `projinfo` dumps, normalized), `run.log`, `BENCH/*`, `*.provenance` and `metadata.json` (judged after
    normalization: the run directory replaced, the TIMINGS, SPLIT TIMINGS and THREAD-SECONDS sections of
    BENCH dropped, and lines with a wall time, RSS, PID or byte count dropped; all NORM_SAME on the base),
    `EXIT_CODE_JOB_1.job` and `time.log` (host, PID, timing only). The base's first run of each case is
    kept as the "before" of later identity checks: `scratch/keep/identity_base/id_base_{cont,prob_neigh}_1`.
  - **Per-client baselines (R1).** A scratch copy of the base with stdout tags after each L-BFGS-B solve
    (evaluations, iterations from `isave(30)`, convergence, wall time of the solve), the Cartesian
    no-reliable-update causes and three accuracy prints (`scratch/instr/patch_instr.py`, tree
    `scratch/instr/rel`, deleted after the runs; rebuilt from `base_src` and the script) changes no number: its identity runs are byte-identical to the base's
    (`scratch/logs/cmp_base_instr_{cont,prob_neigh}.txt`) and its suite logs differ from the base's only in
    PIDs and timings. Counts were identical in four runs of the unit suites
    (`scratch/runs/bl_instr_{1..4}/`); times below are the quiet repeats 3 and 4 (one job of mine at a
    time; load 3.1-6.3 from that job and the OS daemons), 1 thread, so seconds are thread-seconds; the
    `unit_pftc_align2D3D` count includes its thread-consistency test, which sets its own team (capped at 4 by `OMP_THREAD_LIMIT`).
    Summary script `scratch/summarize_baseline.py`.

    | Client (fixture) | Solves | Evals median / mean | Iterations median / mean | Converged | s per solve, median (run 3 / 4) | Accuracy at the base (floor) |
    | --- | --- | --- | --- | --- | --- | --- |
    | `pftc_shsrch_grad` shift, alternating (`unit_pftc_align2D3D`) | 3086 | 11 / 16.60 | 3 / 3.76 | 85.6 % | 6.5e-5 / 6.5e-5 | the shift-prior and in-plane-state floors (0.25, 0.1, 0.15 px) pass |
    | same (`unit_numerics`) | 264 | 12 / 14.21 | 4 / 4.27 | 98.1 % | 3.2e-5 / 3.4e-5 | suites pass |
    | `pftc_shsrch_grad` shift, fixed angle (`unit_pftc_align2D3D`, case e) | 1 | 11 | 4 | yes | 4.1e-5 / 4.3e-5 | 0.126 px rms from the truth shift at the grid angle, 2 deg off (no floor) |
    | `pftc_shsrch_grad` joint (`unit_pftc_align2D3D`) | 6141 | 19 / 24.18 | 6 / 5.87 | 95.2 % | 2.49e-4 / 2.50e-4 | case e: angle error 2.0000 → 0.0001 deg, shift 0.0001 px rms (0.25 px) |
    | `ftexp_shsrch` (`unit_numerics suite=shift_search`) | 100 | 12 / 12.36 | 4 / 4.08 | 100 % | 1.4e-5 / 1.1e-5 | lowest correlation 1.00000 (0.999), largest shift error 1.43e-6 px (0.05) |
    | `ctf_estimate_cost4Dcont`, 3 variables (`mini_stream suite=6vxx`, 8 threads) | 36 | 12 / 13.22 | 6 / 7.72 | 100 % | 2.50e-2 / 2.78e-2 (wall) | defocus error per micrograph 0.0068, 0.0020, 0.0010, 0.0009, 0.0136, 0.0003 µm (0.15) |
    | `quality_logistic`, toy, 4 coefficients (`unit_numerics`) | 3 | 27 / 27 | 8 / 8 | 100 % | 1.69e-4 / 1.68e-4 | cos 0.99702 (0.99), clean rows 0.97368 (0.95), objective 0.27679 |
    | same, `cavg quality nano`, 23-107 coefficients | 7 per size | 35-41 / 34.9-38.9 | 9-10 / 8.9-9.7 | 100 % | 1.8e-4 - 6.9e-4 | suite passes |

    No unit tester reaches `ctf_estimate_cost4Dcont` (`simple_ctf_tester` covers the model only), so its
    fixture is the existing `mini_stream` test (simulated micrographs, distributed CTF fit; decision under
    delegation, R1); both runs gave identical CTF counts and defocus errors
    (`scratch/runs/mini6vxx_instr_{1,2}/`, 42 s each). On the same fixture (refine2D, 8 threads) the polar
    searches ran 1080 shift solves (median 7 evaluations) and 675 joint solves (median 13; mean 16.96 and
    17.06: not reproducible at 8 threads). The `prob_neigh` identity case at one thread ran 25 000 shift
    solves (median 7, mean 6.59 evaluations, 3.6e-4 s) and 2567 joint solves (median 12, mean 16.69,
    2.9e-3 s) (`scratch/runs/id_instr_prob_neigh/SIMPLE_SUBPROC_OUTPUT`).
    **Dose series** (`unit_numerics suite=dose_series`, 295 checks): 1043 Guinier fits; the hand-written LM
    exited at iteration 2 in 649, 3 in 375, 4-33 in 18, and reached the 50-iteration cap once (identical in
    all runs); the suite's numbers: critical-exposure scatter / mean reported error 0.9115, neighbouring
    windows 0.9811, relative B 0.9424; worst relative errors of the model critical exposure 6.141e-4 and
    7.831e-4; exported curve ln c and B 2.119e-4 and 1.037e-2; neighbouring-window exposure 1.502e-2;
    accelerating fading 4.792e-3 and 4.980e-1 (`scratch/logs/dose_series_numbers_base.txt`).
    `dose series manifest` is a `unit_reconstruction` suite and never calls the fit.
  - **Cartesian no-reliable-update outcomes** (tags at the veto `if( .not. identifiable )` and at the
    damped solve's `if( .not. reliable )` of `refine_joint`, and at the two tests of `refine_shift`):
    identity case `refine=cont`, 2500 solves (BENCH: all accepted): 0 veto, 0 failed damped solve;
    `unit_cart_align3D`: 1 veto (`test_invalid_and_unobservable_inputs`, the zero reference that must give
    status 3), 0 failed damped solves; `lib_cart_align3D`, 8000 solves on both routes: 0 and 0, no
    shift-stage failure (statuses: joint all accepted; `shift_then_joint` euclid 997+3 and 999+1 accepted
    plus no improvement). No damped-solve failure on the base: decision 3's stop does not apply.
  - **`cont_route` references** (dead-code probe: the script's `%` form needs a procedure binding, so the
    components were renamed in their module families with its own helpers, `scratch/probe/probe_components.py`;
    four builds, `scratch/logs/probe_build{,2,3,4}.log`, the fourth clean). Families: `simple_parameters`
    (`.f90` 144, 669; `_parse.f90` 225; `_phases.f90` 765-770, 1104) for `cont_route` and
    `l_cont_shift_first`; `simple_cartft_pose_opt` (`.f90` 75, 190-230, 285, 301, 359) for `shift_first`.
    Users outside them: production `strategies/search/simple_strategy3D_cont.f90:143`; tests
    `cftc/simple_cartft_pose_opt_tester.f90:592, 715` and
    `strategies/search/simple_strategy3D_cont_1jyx_tester.f90:254, 255, 258, 259`. The key as a string:
    `params/simple_parameters_tester.f90:86, 91`, `commanders/test/simple_commanders_test_highlevel_cont.f90:717`.
    All are in Phase 1's row; `ui/simple_ui_visibility_tester.f90` names none of them.
  - **Findings for later phases.** L-BFGS-B leaves `ospec%niter` at 0 (only simplex sets it) and counts
    `nevals == ngevals`; Phase 3's tables compare evaluations, and an L-BFGS-B iteration count would need
    `simple_opt_lbfgsb.f90`, which is not in the table. `ctf_estimate_cost4Dcont` has no unit test; Phase 3
    adds its oracle and comparison in `simple_ctf_tester.f90`.
  - Exit: recorded here. What it means, in short: the base builds and gates clean, both identity cases
    are deterministic at one thread with the exclusions above, the Cartesian route has never needed the
    damped-solve failure path on these fixtures, and every client has counted, repeatable L-BFGS-B
    baselines to compare against.
- **2026-11-17, Phase 1 (identical, the shift-only path aside), arun `lm_optimizer`.** FWM-S116118, toolchain as
  in Phase 0; base `aae467bc7`. Paths below are under `/Users/elmlundho/agent_work/arun/lm_optimizer`.
  Files changed: `src/utils/math/simple_lm_solver.f90` (new: `lm_problem`, `lm_options`, `lm_outcome`, `lm_solve`,
  `lm_scale_system`, the `LM_*` outcomes), `src/utils/math/simple_linalg.f90` and its tester (`spd_solve_rel`,
  `solve_pose_cholesky` generalised to n with its expressions; `test_spd_solve_rel`, 6 checks),
  `cftc/simple_cartft_pose_opt.f90`, `_system.f90`, `_tester.f90` (`cartft_pose_problem`, `cartft_shift_problem`,
  `refine_joint` and `refine_shift` on `lm_solve`, `refine_shift` public, `CARTFT_*` defined as the `LM_*` values,
  `pose_data_observable` the veto; `build_pose_lm_system`, `prepare_pose_data_system`, `scale_pose_system`,
  `apply_pose_parameter_mask`, `solve_pose_cholesky`, `increase_lm_damping`, `rejection_terminal_status` and the LM
  constants gone from `cftc/`), the `cont_route` retirement (`params/simple_parameters{,_parse,_phases,_tester}.f90`,
  `strategies/search/simple_strategy3D_cont.f90`, `simple_strategy3D_cont_1jyx_tester.f90`,
  `commanders/test/simple_commanders_test_highlevel_cont.f90`), `doc/policies/3D/refine3D_policy.md`, skills
  `simple-main-cftc` and `simple-refine3d`, this note. All are in Phase 1's rows; `ui/simple_ui_visibility_tester.f90`
  was not needed (it names only the older `pose_cont_route`).
  - **Builds and gates.** Release and Debug, no compiler warning; fast gate 15/15 on both (Release 16.4 s, Debug
    55.9 s: `scratch/logs/fastgate_p1x_release_final.log`, `fastgate_p1x_debug.log`; the gate's one warning names the
    new, untracked `simple_lm_solver.f90`, as nothing is committed). Identity trees are built from an export of the
    working tree without the repository metadata (`scratch/build_export.sh`, trees `scratch/p1x_{release,debug}`), as
    the base was: a build from the checkout writes the commit into `run.log` (`SIMPLE Git Commit:`), which the base,
    built from `git archive`, leaves empty (seen on a first, checkout-built tree: `scratch/logs/cmp_p1rel_cont.txt`).
    Phase 0's exclusions are unchanged; `scratch/compare_runs.sh` now also maps the original run directory of a base
    run moved to `scratch/keep/` to `RUNDIR` (the declared run-directory exclusion).
  - **Exit 1, the Cartesian identity: met.** Against the Phase 0 base runs, at one thread, on Phase 0's identity
    outputs: Release `refine=cont` 63 files byte-identical, the 12 excluded files identical after normalization
    (`scratch/logs/cmp_p1x_rel_cont.txt`); Debug `refine=cont` the same, against a Debug base pair run here and found
    deterministic (63 + 12, `cmp_basedbg_cont.txt`; Phase 1: `cmp_p1x_dbg_cont.txt`; the base Debug run is kept in
    `scratch/keep/identity_base/id_basedbg_cont_1`); Release `refine=prob_neigh` the same (`cmp_p1x_rel_prob_neigh.txt`).
  - **Exit 2, the cart suites: met.** `unit_cart_align3D` and `lib_cart_align3D` print the base's numbers line for
    line, less the retired route's lines (`scratch/runs/suites_p1_rel/` against `scratch/runs/suites_rel1/`).
    Assertions 1183 → 1176 (the shift-then-joint half of E16 and the third case of the committed-residual test removed;
    E16's pose-update assertion moved to its joint case) and 106 → 52 (the `shift_then_joint` route, 52 checks, and
    the joint route's two checks of the `cont_route` default). `test_prior_identifiability` passes with its assertions
    unchanged, its calls renamed (`lm_scale_system` and `pose_data_observable` for `prepare_pose_data_system`,
    `spd_solve_rel` for `solve_pose_cholesky`); `test_invalid_and_unobservable_inputs` passes unchanged (the zero
    reference still ends in `CARTFT_NO_RELIABLE_UPDATE` by the veto).
  - **Exit 3, the shift-only cases: met.** `test_shift_solver` passes at its tolerances (recovery within `SHIFT_TOL`
    2e-3 px, fewer than 20 iterations, step bound kept; an exact shift reports no improvement with no attempt).
  - **Exit 4, the searches: met.** `rg` finds none of the moved helpers or LM constants in `cftc/`, and `rg -w` finds no
    `cont_route`, `l_cont_shift_first` or `shift_first` in `src/` (a substring search also hits the retired key
    `pose_cont_route` in `ui/simple_ui_visibility_tester.f90`, which checks that no program offers it).
  - **Other runs.** `unit_numerics` 2066 assertions (2060 + `test_spd_solve_rel`), `unit_image`, `unit_pftc_align2D3D`
    and `fills` (48 rows, no differing entry) pass. `test=cont_refine3D_1jxy nthr=8` (Release; its commander changed):
    116 PASS, 0 FAIL, 277 s, load 11.6-13.0 from my other runs (`scratch/runs/hl_cont_p1x/`); stages a-c, now on the
    joint route, at E25's floors: median rotation 0.958 / 0.229 / 0.958 deg, shift 0.067 / 0.058 / 0.067 px, basin
    1.0, c equal to a.
  - **Decisions under delegation.**
    1. The engine has one more hook than decision 1 lists, `negligible_step` (default: never). `refine_joint` also
       stops after an accepted step whose rotation and shift norms are both below 1e-8; a generic step norm cannot
       reproduce a per-group test, and decision 3 extracts the step control verbatim. The Cartesian problems override it.
    2. `refine_joint`'s seed handling stays in the Cartesian wrapper (the transaction always passes seed terms); the
       engine evaluates the start itself only without them (`eval` evaluates the trial point, which is the start
       until the first `propose`).
    3. `refine_shift` keeps its residual-floor early exit (an exact shift: `CARTFT_NO_IMPROVEMENT`), and its
       `propose` bounds trials to `trs` from the input shift (decision 10); it has no observability veto (decision
       6: the generic engine vetoes nothing).
    4. N20's stages i and j (a and b on the joint route) were deleted: with the seam retired they are the same
       commands as a and b.
  - **Finding for Phases 2 and 3 (decision 6 does not hold as written).** With `d_i = max(C_ii, sqrt(eps)*max|C|,
    eps**2)` and `spd_solve_rel`'s pivot floor `sqrt(eps)*max|C + mu*diag(d)|`, a coordinate whose curvature is
    below `sqrt(eps)*max|C|` (a null direction: `C_ii = 0`) never passes the pivot test, at any μ: its damped pivot is
    `mu*sqrt(eps)*max|C|`, the floor about `(1+mu)*sqrt(eps)*max|C|` (checked numerically for μ from 1e-3 to 1/eps).
    The damped solve then fails every iteration and the engine ends with no improvement after 8 rejections, instead
    of giving that coordinate a zero step and moving the others. The Cartesian route never reaches it (the veto
    rejects such a data system first, and Phase 0 counted no damped-solve failure), so Phase 1's identity is
    unaffected; Phase 3's isotropic CTF start (zero angle Jacobian) would. A remedy that keeps the Cartesian identity:
    after the observability test, give a coordinate with zero scaled gradient and zero scaled curvature the inactive
    mask's unit diagonal for the damped solve. It is not applied here: it is a change to decision 6, for the phase
    whose test exercises it (Phase 2's L5, whose curvature diagonal is nonzero, does not).
  - What it means, in short: the library has one LM engine, the Cartesian route runs on it byte for byte in Release
    and Debug, `shift_then_joint` is gone, and the engine's null-direction handling needs a decision before a client
    with an exactly singular curvature (the isotropic CTF fit) uses it.
- **2026-11-17, Phase 2 (tests-only), arun `lm_optimizer`: STOPPED at L2 and L3.** FWM-S116118, toolchain as in
  Phase 0; base `aae467bc7`. Paths below are under `/Users/elmlundho/agent_work/arun/lm_optimizer`.
  Files changed: `opt/simple_opt_lm.f90` (new: `opt_lm`, the private `spec_problem` with `eval`, `free`,
  `limit_step`, `propose`, `accept`), `opt/simple_opt_spec.f90` (`fdfcurvfun_8`, `set_fdfcurvfun_8`,
  `eval_fdfcurv_8`, `'lm'` in `specify`), `opt/simple_opt_factory.f90` (`'lm'`), `opt/simple_opt_tester.f90` (L1-L5,
  the plain-clipping control), the `simple-main-opt` skill, this note. Release build clean, no warning
  (`scratch/logs/build_p2_release.log`).
  - **Results** (`unit_numerics suite=optimisers`, `scratch/logs/p2_optimisers_suite.log`; 95 checks, 3 fail):
    L1 passes (error below 5e-9; 30 evaluations against L-BFGS-B's 100 on the same cost). L4 passes (4 callbacks,
    each followed by one counted evaluation at its point; the solve continues on the changed objective). L5 passes
    (f < 1e-12; x - y keeps its value within 1e-12). L2's plain-clipping control passes (it stalls at the start,
    status bound-rejected). **L2 fails:** one accepted step to (0, -0.01998), then the next proposal predicts less
    than 1e-6·|f| (f = 2.304 at the minimum) and the solve ends unevaluated: error 2.0e-5 (criterion 1e-6),
    projected gradient 1.0e-4 (criterion 1e-8). **L3 fails its count:** no step accepted and the start returned,
    but after 7 rejections; the 8th proposal predicts less than 1e-6·|f| and ends the solve unevaluated
    (8 evaluations, not 9). A double-precision emulation of the engine predicted all of this before the build
    (`scratch/p2/emulate.py`).
  - **Mechanism.** Decision 3's stops relative to |f| (predicted reduction below 1e-6·|f|, accepted relative
    reduction below 1e-6) end a solve when the remaining gain is small against the cost itself. L2's constrained
    minimum has a large cost, so they fire about 2e-5 from it, long before the projected-gradient test that
    decision 4 counts on. In L3 the damping grows by 4, 8, 16, ... per rejection, and after 7 rejections a
    proposal's predicted reduction is below 1e-6 of a least-squares cost. With the two tolerances as a
    parameter (`scratch/logs/p2_tolerance_variants.txt`): L3 reaches 8 rejections from 1e-9 down; L2's 1e-6
    accuracy needs 1e-12 or less, and its projected gradient below 1e-8 holds only with no relative stop (0:
    error 2.5e-12, projected gradient 1.3e-11); L1 is unchanged by any of them.
  - **Stop.** Plan section 4: "An LM test passes only with a constant other than `refine_joint`'s." Not done:
    the mutation tests of L2's control and L3 (they follow the ruling), the fast gate (it fails on this suite),
    Debug. `rg "'lm'"` in `src/` finds no caller outside `opt/` (its definition and the factory entry).
  - **For Hans to decide** (one ruling covers both): (a) L2 and L3 restated at what `refine_joint`'s constants
    give (L2: error within the predicted-reduction stop's reach, about sqrt(2e-6·|f|/λ), and no projected-gradient
    floor; L3: no step, the start returned, no count); or (b) the engine's two relative tolerances become options,
    `refine_joint`'s values by default (the Cartesian identity is unchanged), and the factory passes smaller ones
    (0 for L2's criteria, or a value from `opt_spec`), which then also decides what `factr`/`pgtol` mean under `lm`.
    The Phase 1 finding on decision 6 (a null coordinate never passes the damped solve's pivot floor) needs a
    ruling before Phase 3's isotropic CTF oracle as well.
- **2026-11-18, Phase 2 completed under R3 and R4, arun `lm_optimizer`** (supersedes the stop above). FWM-S116118,
  toolchain as in Phase 0; base `aae467bc7`. Paths below are under `/Users/elmlundho/agent_work/arun/lm_optimizer`.
  Files changed since the stop: `src/utils/math/simple_lm_solver.f90` (R3: `lm_options%predicted_tol`,
  `accepted_tol`, `stationary_tol`, defaulting to 1e-6, 1e-6 and 1e-8; R4: after the observability test and
  the stationarity test, a coordinate with exactly zero scaled gradient and curvature diagonal gets the
  inactive mask's unit diagonal), `opt/simple_opt_lm.f90` (both relative stops at `factr`·ε, the gradient stop
  at `pgtol` on the free gradient), `opt/simple_opt_tester.f90` (tolerances per test, L6, the control at L2's
  tolerances), the `simple-main-opt` skill, this note (R3 and R4 in section 3, decisions 3 and 6, section 4,
  the Phase 2 row; Phase 2 added to `simple_lm_solver.f90`'s row, as R4's "Where R3 and R4 go" says).
  - **Exit 1, L1-L6: met** (`unit_numerics suite=optimisers`, 99 checks: `scratch/runs/p2_opt/optimisers_2.log`
    Release, `scratch/logs/p2_optimisers_debug.log` Debug). The tolerances each test sets through `specify`
    (R3; chosen with the engine's emulation `scratch/p2/ruled.py` before the build): L1, L4, L5, L6 `factr` 1e7,
    `pgtol` 1e-8; L2 `factr` 0 (no relative stop, L-BFGS-B's meaning), `pgtol` 1e-10; L3 `factr` 1e7, `pgtol`
    1e-5 (L-BFGS-B's defaults). L1: error below 5e-9, 30 evaluations (L-BFGS-B 100 on the same cost, reported).
    L2: (0, -0.02) within 1e-6 (the y error rounds to 0 at the printed 7 digits), projected gradient 1.2e-11,
    4 evaluations; its control (the engine with the box clip and no free hook, at L2's tolerances) does not
    reach it (y = -0.0077 after 65 accepted steps of 100 iterations). L3: the start returned bit for bit after 8
    rejections (9 evaluations). L4: 4 callbacks, each followed by one counted evaluation at its point, and the
    solve finishes on the changed objective. L5: f < 1e-12, x - y within 1e-12 of its start. L6: f < 1e-12,
    y exactly its start value, the engine's status accepted.
  - **Mutation tests** (`scratch/p2/plant.py`, `run_plant.sh`, on the export copy only; logs
    `scratch/logs/mut2_*.log`): P1, `opt_lm` without its free hook (plain clipping): L2's two checks fail;
    P2, the engine accepts every evaluated trial: L3's three checks fail (and L1); P3, the control's solve given
    the free rule: the control fails. Restored: 99 of 99. A first round ran with the plants accumulating (rsync
    kept the restored files' older time stamps, so make did not rebuild them) and showed that the control, then
    at the engine's default tolerances, could not fail (a solve with the free hook also stopped 2e-5 short); the
    control now runs at L2's tolerances, and the harness touches the restored files.
  - **Exit 2, no production caller: met.** `rg "'lm'"` in `src/` outside the testers finds only the optimiser's
    definition, the spec and the factory entry.
  - **Cartesian identity rerun (R4's condition): met**, against the Phase 0 base, at one thread: Release
    `refine=cont` 63 files byte-identical + 12 normalized (`scratch/logs/cmp_p2_rel_cont.txt`), Release
    `refine=prob_neigh` the same (`cmp_p2_rel_prob_neigh.txt`), Debug `refine=cont` the same against the base
    Debug run (`cmp_p2_dbg_cont.txt`). `unit_cart_align3D` and `lib_cart_align3D` print Phase 1's lines
    exactly (`scratch/runs/suites_p2_rel/` against `suites_p1_rel/`).
  - **Gates.** Release and Debug builds without warnings (`scratch/logs/build_p2r_release.log`,
    `build_p2r_debug_final.log`); fast gate 15/15 on both (`fastgate_p2_release_final.log`,
    `fastgate_p2_debug_final.log`); `unit_numerics` 2087 assertions, `unit_image`, `unit_pftc_align2D3D` and
    `fills` (48 rows, no differing entry) pass.
  - **Decisions under delegation.** L1, L4, L5 and L6 also set `factr`/`pgtol` explicitly (R3 names L2 and
    L3): at the spec's default `pgtol` 1e-5 the emulation stops L1 1.5e-5 from (1, 1), outside its 1e-6; the
    criteria are unchanged. L6 checks the engine's status with a tester problem on `lm_solve` (the factory
    exposes none). R4's mask is applied after the stationarity test, so it changes only the damped system.
  - What it means, in short: the factory has a Levenberg-Marquardt optimiser that runs on the library's one
    engine under each client's own `factr`/`pgtol`, handles bounds by freezing coordinates and null directions
    by a zero step, and is caught by its tests when either is broken; the Cartesian route is still byte-identical.
- **2026-11-18, Phase 3 (tests-only for the factory clients, production identical; result-changing for the dose
  series), arun `lm_optimizer`.** FWM-S116118, toolchain as in Phase 0; base `aae467bc7`. Paths below are under
  `/Users/elmlundho/agent_work/arun/lm_optimizer`. Files changed: `pftc/simple_polarft_calc.f90` and
  `simple_polarft_corr.f90` (`gen_shift_curv_for_rot_8_aln`, `gen_joint_curv_at_angle`, `ring_partner`),
  `pftc/simple_pftc_shsrch_grad.f90` (`SHSRCH_SOLVER`, `solver` in `new_alternating`/`new_fixed`/`new_joint`, the
  `'lm'` callback, `lm_cost_to_score_cost`, the test seam's interface) with its new tester submodule
  `simple_pftc_shsrch_grad_tester.f90` (the seam `shsrch_lm_terms`), `pftc/simple_pftc_inplane_tester.f90`;
  `image/simple_ftexp_shsrch.f90` and its tester, `image/simple_ft_expanded.f90` (`corr_shifted_curv_8`);
  `ctf/simple_ctf_estimate_cost.f90` (`CTF_CONT_SOLVER`, `solver` in `init`, the curvature in `fdf`) with its new
  tester submodule `simple_ctf_estimate_cost_tester.f90`, `ctf/simple_ctf_tester.f90` (calls it);
  `quality/simple_quality_logistic.f90` and its tester; `src/utils/filter/simple_dose_series.f90` (the Guinier fit on
  `lm_solve`); this note. The two new testers are covered by the table's rule that a listed `x.f90` covers
  `x_tester.f90`; `simple_ft_expanded.f90` is its "if needed" row.
  - **Exit 1, the curvature oracles: met, no client excluded.** At a zero residual the curvature equals the
    central-difference Hessian of the analytic gradient and is symmetric positive semi-definite. Tolerances fixed in
    the testers before their runs: `quality_logistic` 1e-6 relative (guessed; the Fisher matrix plus the ridge is the
    exact Hessian, checked at the start and at a point away from it): 5.0e-11 and 2.4e-12; `ftexp_shsrch` 1e-3
    (the tester's existing gradient check): 1.2e-7 (peak of an integer-shifted square); `ctf_estimate_cost4Dcont`
    1e-3 (guessed), data proportional to the model, 3 and 4 parameters: 5.7e-8 and 7.4e-8; `pftc_shsrch_grad` at the
    gradient check's `GRAD_FD_RTOL` 1e-2, the particle the reference rotated by a whole grid index and shifted:
    shift mode euclid 2.1e-8 and cc 2.1e-8, joint mode euclid with the prior centred at the solution (nonzero centre,
    angle active, so the 2w[I,-c']^T[I,-c'] term is in the matrix) 7.7e-6, joint cc 8.7e-5. CTF from an isotropic
    start (dfx = dfy = 2.02 µm against 2.00): the defocus improves (to 1.9946, dfx = dfy throughout) and the angle
    stays within 2e-6 degrees (tolerance 1e-4, guessed); the null-coordinate rule R4 is what lets it move.
  - **Exit 2, one comparison table per client: met.** Release, one thread, so seconds are thread-seconds; counts
    identical in two runs (`scratch/runs/cmp_{1,2}/`), times of run 1 / run 2; load 5.5-8 during both (the session's
    earlier jobs fading and desktop activity on the Mac): times are secondary and Phase 4 re-times them quietly.
    The L-BFGS-B rows run the production code unchanged (identical counts to Phase 0 where the fixture is
    Phase 0's: `ftexp_shsrch` 12.36, logistic 27). L-BFGS-B's iteration count is not reported (it leaves `niter` at
    0, Phase 0 finding); `'lm'`'s is.

    | Client, fixture | Solver | Evaluations per solve | Iterations | s per solve | Accuracy |
    | --- | --- | --- | --- | --- | --- |
    | `quality_logistic`, toy, 4 coefficients | lbfgsb | 27 | (8 in Phase 0) | 1.8e-4 / 1.9e-4 | cos 0.997017, clean 0.973684, obj 0.27678898 |
    | | lm | 7 | 7 | 1.2e-4 / 9.6e-5 | cos 0.997016, clean 0.973684, obj 0.27678898 |
    | `ftexp_shsrch`, 100 sub-pixel shifts | lbfgsb | 12.36 | (4.08 in Phase 0) | 1.25e-5 / 1.16e-5 | corr 1.00000, max shift error 1.4e-6 px (floor 0.05) |
    | | lm | 3.64 | 3.64 | 7.6e-6 / 6.9e-6 | corr 1.00000, max shift error 3.3e-4 px |
    | `ctf_estimate_cost4Dcont`, 8 starts inside the basin, 3 parameters | lbfgsb | 57.38 | | 1.06e-2 / 1.07e-2 | defocus error mean 0.00184, max 0.00186 µm; angle 0.037 deg |
    | | lm | 29.88 | 29.50 | 6.1e-3 / 5.8e-3 | defocus error mean 0.00183, max 0.00201 µm; angle 0.040 deg |
    | `pftc` euclid, alternating (3 starts) | lbfgsb | 43.67 | | 1.7e-4 / 1.9e-4 | shift error mean 0.089, max 0.178 px; angle 0.5 deg |
    | | lm | 8.33 | | 6.3e-5 / 7.2e-5 | shift error mean 0.045, max 0.045 px; angle 0.5 deg |
    | `pftc` euclid, fixed angle | lbfgsb | 12.33 | | 3.5e-5 / 4.3e-5 | shift error 0.178 px (grid angle 2 deg off) |
    | | lm | 4.33 | | 1.7e-5 / 1.9e-5 | shift error 0.178 px |
    | `pftc` euclid, joint | lbfgsb | 19.00 | | 3.3e-4 / 3.8e-4 | shift error max 0.0003 px; angle 0.001 deg |
    | | lm | 4.00 | | 2.2e-4 / 2.3e-4 | shift error max 0.0002 px; angle 0.001 deg |
    | `pftc` cc, alternating | lbfgsb | 20.00 | | 1.0e-4 / 1.2e-4 | shift error mean 0.089, max 0.178 px |
    | | lm | 9.00 | | 6.4e-5 / 7.5e-5 | shift error mean 0.045, max 0.045 px |
    | `pftc` cc, fixed angle | lbfgsb | 10.33 | | 3.7e-5 / 3.5e-5 | shift error 0.403 px (grid angle 4.5 deg off) |
    | | lm | 4.67 | | 2.0e-5 / 2.2e-5 | shift error 0.403 px |
    | `pftc` cc, joint | lbfgsb | 19.00 | | 3.6e-4 / 3.8e-4 | shift error 0.0002 px; angle 0.000 deg |
    | | lm | 4.00 | | 2.1e-4 / 2.4e-4 | shift error 0.0002 px; angle 0.000 deg |

    Phase 0's spread on the same machine: counts identical in four runs, times within about 20% (its table above).
    A `pftc` "solve" is one `minimize` call (the alternating search runs several L-BFGS-B or LM solves inside it).
    Decision 8 is not applied here (Phase 4).
  - **Exit 3, the polar identity: met** (`scratch/logs/cmp_p3_rel_prob_neigh.txt`: Release `refine=prob_neigh`
    against the Phase 0 base, 63 files byte-identical, 12 identical after normalization); the Cartesian identity also
    still holds (`cmp_p3_rel_cont.txt`, 63 + 12). `polar_matching_bench suite=fills nthr=8`: 48 rows, no differing
    entry (`scratch/runs/p3_fills/polar_fill_identity.tsv`). The CTF fit's production numbers are unchanged:
    `mini_stream suite=6vxx` gives Phase 0's six defocus errors digit for digit (`scratch/runs/mini6vxx_p3/run.log`).
  - **Exit 4, the dose series at its floors: met.** The Guinier fit runs on `lm_solve` (2 parameters, unbounded,
    relative stops 1e-10 as the old loop's, 50 iterations); `dose series` passes its 295 checks and prints Phase 0's
    numbers to their printed precision (`scratch/logs/dose_series_numbers_p3.txt` against `_base.txt`), the same
    windows are reported unfitted; `dose series manifest` (in `unit_reconstruction`, which does not run the fit)
    passes its 87.
  - **Gates.** Release and Debug without warnings (`scratch/logs/build_p3_release3.log`, `build_p3_debug2.log`); fast
    gate 15/15 on both (`fastgate_p3_release.log` 16.4 s, `fastgate_p3_debug.log` 50.7 s); `unit_numerics` 2102,
    `unit_image` 982, `unit_pftc_align2D3D` 6335, `unit_reconstruction` 1378 assertions pass in Release and Debug
    (`scratch/runs/p3b/`, `p3_dbg/`; the "Error termination" lines in the Debug logs are the suites' intended child
    failures, the same counts as the base Debug build's, `scratch/runs/base_dbg_num/`). The fast gate's
    `check_tester_linkage` first refused the seam calling a private procedure; it now calls through the spec's
    pointer. Production: each client's solver constant is `'lbfgsb'`, no production call passes a solver, `'lm'` has
    no production caller.
  - **Decisions under delegation.**
    1. Curvatures. `ftexp_shsrch`: (1e8/den)·Σ w|T| t t^T (T = A·conj(B)), the Hessian at a zero residual and the
       Gauss-Newton matrix of the equivalent least squares, constant over the shift: computed once per `'lm'` solve.
       CTF: the normalised model's Gauss-Newton matrix (J^T J)/P - (J^T a)(J^T a)^T/P^2 plus the astigmatism
       penalty's exact Hessian, without `cartft_calc`'s `curvature_scale` slope (a Cartesian heuristic for noisy data;
       it is 1 at a zero residual). Polar shift mode: euclid (2/W) Σ w_k a a^T |R CTF|^2, cc the same with k weights
       over sqrt(D C); joint mode: the three derivative sections rotated to the fractional index by an inverse FFT
       done with the existing forward plan (conjugate in, conjugate out), summed with CTF^2 over the full ring;
       cc normalised with the denominator series. Logistic: the Fisher matrix plus the ridge.
    2. The polar shift mode under `'lm'` minimises the euclid loss L/W + w|xy - c|^2 (the log of the score with its
       prior, decision 6's premise) and converts the result back to the L-BFGS-B cost -exp(-loss) where `minimize`
       compares costs; a cc shift prior under `'lm'` stops (`THROW_HARD`): there is no least-squares form of cc·g,
       and no production route combines them with `'lm'`.
    3. The production `fdf` of the CTF cost gains an optional curvature output in its loop; `mini_stream` shows its
       L-BFGS-B numbers unchanged. The polar and `ftexp` curvatures are separate routines (production paths untouched).
    4. Oracles of private procedures run in tester submodules: the CTF cost's new tester submodule, and for the shift
       search a seam (`shsrch_lm_terms`, interface in the parent, body in its tester submodule, so it exists only in
       test builds) that the in-plane tester calls.
    5. Fixtures: the CTF tests start inside the cost's capture basin (0.02 µm and 5 degrees at most: at 5 A the phase
       moves about 25 rad per µm, so the first run's 0.15-0.25 µm starts ended in neighbouring local minima under both
       solvers); `ftexp_shsrch`'s comparison replays the same 100 shifts for both solvers.
    6. Exact symmetry: the logistic and `ftexp` curvatures mirror their upper triangle ((w x_i) x_j and (w x_j) x_i
       round differently; the first run failed the exact-symmetry check).
  - What it means, in short: every factory client can run under `'lm'` with a curvature verified against its
    Hessian; on these fixtures `'lm'` needs 2-7 times fewer evaluations than L-BFGS-B at comparable accuracy (better
    in the polar alternating search, slightly worse at the tail in the CTF fit, 0.00201 against 0.00186 µm, and in
    `ftexp_shsrch`, 3e-4 px against 1e-6, far inside its 0.05 floor); production is byte-identical; and the
    dose-series fit runs on the library's engine. Phase 4 judges these by decision 8 on quiet timings and real data.
