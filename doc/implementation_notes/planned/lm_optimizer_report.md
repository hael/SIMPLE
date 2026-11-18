# One Levenberg-Marquardt for the library: interim report (Phases 0-3)

Status, 2026-11-18: interim. The run `lm_optimizer` (arun, FWM-S116118, unattended) carried out Phases 0-3 of
`lm_optimizer.md` on base `aae467bc7` (ruling R2). Phase 4 (decision 8's verdicts, default switches, real-data A/B)
is a later run with Hans and extends this report. Nothing is committed; the change is the working tree. The plan's
Progress section holds the detail; paths of evidence are under `/Users/elmlundho/agent_work/arun/lm_optimizer`.

## 1. What was learned

- **One engine carries every LM client, the Cartesian route byte for byte.** `refine_joint` extracted into
  `lm_solve` reproduces the base's `refine=cont` outputs in Release and Debug (63 files byte-identical, the declared
  exclusions identical after normalization), rechecked after every later engine change (Phases 1, 2, 3).
- **The step control as extracted did not suit every client** (Phase 2 stop, rulings R3 and R4):
  - the stops relative to |f| (1e-6) end a solve far from a minimum whose cost is large (L2 stopped 2e-5 from it);
    they are now options, `refine_joint`'s values by default, and the factory maps `factr`/`pgtol` onto them;
  - a coordinate with exactly zero curvature never passed the damped solve's relative pivot floor at any damping
    (the Phase 1 finding); such a coordinate now takes a zero step (R4). A coordinate with tiny but nonzero
    curvature (below sqrt(eps) of the largest) still cannot: see 6.
- **Every factory client's Gauss-Newton curvature equals its Hessian at a zero residual** (oracles: 2e-12 to 9e-5
  relative), including the polar joint mode with the angle active and a nonzero prior centre.
- **On the testers' fixtures `'lm'` needs 2-7 times fewer evaluations than L-BFGS-B** at comparable accuracy
  (Phase 3 table): logistic 7 against 27, `ftexp_shsrch` 3.6 against 12.4, CTF 30 against 57, polar joint 4
  against 19, polar alternating 8-9 against 20-44. Times (thread-seconds) fall less than the counts because an
  `'lm'` evaluation also forms the curvature (the polar joint one takes two extra FFT batches): 0.37-0.67 of
  L-BFGS-B's per solve across the clients, measured under a load of 5.5-8.
- **L-BFGS-B's iteration count is not observable through the factory** (`niter` stays 0); the comparison uses
  evaluations.
- **The CTF cost is several |sin| periods per 0.1 µm wide at 5 Å** (about 25 rad per µm): the continuous fit is a
  local polish of the grid search; fixtures must start inside its basin.

## 2. What changed, phase by phase

| Phase | Diff | Files | What |
| --- | --- | --- | --- |
| 0 | `review/phase_0.diff` | 1 (this plan) | Measurement: base builds and gates, identity cases and exclusions, per-client baselines, Cartesian outcome counts, the `cont_route` probe; R1, R2 recorded |
| 1 | `review/phase_1.diff` | 17 | `simple_lm_solver` (new), `spd_solve_rel`, `refine_joint`/`refine_shift` on the engine (`cartft_pose_problem`, `cartft_shift_problem`), `cont_route`/`shift_then_joint` retired, policy and two skills |
| 2 | `review/phase_2.diff` | 7 | `opt_lm` (new), `fdfcurvfun_8`, the factory's `'lm'`, L1-L6, R3 tolerances and R4 null coordinates in the engine, `simple-main-opt` skill |
| 3 | `review/phase_3.diff` | 16 (2 new testers, the plan, this report) | Curvatures and solver names of the five factory clients with their oracles and comparisons, the Guinier fit on the engine, this report |

The driver writes `review/phase_3.diff` when this phase ends; `review/full_after_N.diff` holds the cumulative diff.

## 3. Decisions taken under delegation

Each is recorded, with its reason, in the plan's Progress entry of its phase.

- **Phase 0:** the CTF baseline on the existing `mini_stream` fixture (R1: no unit test reaches the continuous CTF
  fit); the identity fixture is the test's seeded project captured before its first write; `nparts=1` in the identity
  commands (the UI requires it); the declared exclusions; the probe of data components through the probe script's
  own helpers (its `%` form needs a procedure binding).
- **Phase 1:** one engine hook beyond the plan, `negligible_step`, for `refine_joint`'s per-group step stop; seed
  handling in the Cartesian wrapper; `refine_shift` keeps its residual-floor exit, bounds trials to `trs`, and has no
  veto; N20's stages i and j deleted (identical to a and b after the retirement); identity trees built from an export
  without repository metadata (the base had no commit hash in its logs).
- **Phase 2:** under R3, L1 and L4-L6 also set `factr`/`pgtol` (their criteria unchanged); L6's status checked on the
  engine directly; L2's plain-clipping control runs at L2's tolerances (at the defaults it could not fail, which the
  first mutation round showed); R4's mask after the stationarity test.
- **Phase 3:** the curvature of each client (plain normalised-model matrix for the CTF, without the Cartesian slope
  scaling); the polar shift mode's `'lm'` objective in log form with its result converted back, and no cc shift prior
  under `'lm'`; tester submodules and a test-only seam for oracles of private procedures; CTF fixtures inside the
  basin; exact symmetry by mirroring.

## 4. Deviations from the plan

- **Phase 2 stopped once** on L2 and L3 (a plan stop: they passed only with constants other than `refine_joint`'s);
  Hans ruled R3 and R4 the same day, which amended decisions 3 and 6 and added L6.
- **Phase 1's exit "`rg` finds no `cont_route` in `src/`"** holds as a word search; a substring search also hits the
  older, already retired `pose_cont_route` in a UI tester that checks no program offers it.
- **N20's stages i and j were deleted** (`simple_commanders_test_highlevel_cont.f90`); `cont_refine3D_1jxy` passed
  116/116 after it.
- **`dose series manifest` is a `unit_reconstruction` suite and never runs the Guinier fit**; it was run, and passes,
  but it does not judge decision 11.
- **The comparison fixtures are the testers'** (R1): the CTF comparison uses the new tester fixture, not Phase 0's
  `mini_stream`; the timings were taken under a load of 5.5-8, partly desktop activity on the Mac.
- **Files outside the table:** none. Two new tester submodules (`simple_pftc_shsrch_grad_tester.f90`,
  `simple_ctf_estimate_cost_tester.f90`) are covered by the table's `x.f90` covers `x_tester.f90` rule;
  `simple_ft_expanded.f90` is its "if needed" row; Phase 2 was added to `simple_lm_solver.f90`'s row by R4.

## 5. Evidence against the baseline

- **Builds and gates:** Release and Debug without warnings in every phase; fast gate 15/15 on both after every phase
  (Phase 3: `scratch/logs/fastgate_p3_{release,debug}.log`).
- **Cartesian identity** (`refine=cont`, Release and Debug) and **polar identity** (`refine=prob_neigh`, Release):
  byte-identical on the declared outputs (`scratch/logs/cmp_p3_rel_{cont,prob_neigh}.txt`; Phase 1 and 2 runs in
  their entries). `polar_matching_bench suite=fills`: 48 rows, no differing entry, every phase.
- **Cart suites** print the base's numbers less the retired route (Phase 1; unchanged in Phase 2).
- **Optimisers suite** L1-L6 pass in Release and Debug (99 checks); the mutation plants on L2, its control and L3
  are each caught (`scratch/logs/mut2_*.log`).
- **Oracles** (Phase 3) and the **comparison table per client** (Phase 3 Progress entry).
- **Dose series:** 295 checks pass and the suite prints the base's numbers to its printed precision.
- **Production unchanged at the end of Phase 3:** each client's solver constant is `'lbfgsb'`, no production call
  passes a solver name, `'lm'` has no production caller; `mini_stream` reproduces Phase 0's defocus errors.

## 6. Left open for the maintainer

- **Phase 4 (R2):** decision 8 per client, with quiet re-timings (these were taken under load) and the real-data
  A/Bs (PorF `local4_best` for the shift searches; the motion-correction and CTF sets Hans names, R1).
- **Reported, not gated:**
  - `ftexp_shsrch` under `'lm'` stops at 3e-4 px against L-BFGS-B's 1e-6 at the same `factr`/`pgtol`;
  - the CTF fit's tail is slightly worse under `'lm'` (0.00201 against 0.00186 µm) on the 8-start fixture;
  - the polar joint curvature recomputes the sections and costs two FFT batches per evaluation; a memo of the
    shift-independent blocks per reference would cut it if `'lm'` wins there.
- **Known limits of the engine:**
  - a coordinate whose curvature is tiny but nonzero relative to the largest still fails the damped solve (R4
    covers exact zeros); the coordinates are unscaled in the factory;
  - `'lm'` with a cc shift prior in the polar shift mode stops by design.
- **Follow-up:** whether the near-null case needs a remedy belongs with Phase 4's real-data reads.

## 7. Commit message (for the maintainer, if the change is taken as one)

```
lm_optimizer, Phases 0-3: one Levenberg-Marquardt engine for the library

simple_lm_solver (lm_problem, lm_solve) owns the LM step control once: the Cartesian pose solve
(refine_joint, refine_shift) runs on it byte-identically, the factory gains 'lm' (opt_lm, fdfcurvfun_8,
L-BFGS-B's factr/pgtol as its stops), the five factory clients gain verified Gauss-Newton curvatures and
a tester-only solver switch, and the dose-series Guinier fit runs on it. cont_route/shift_then_joint
are retired. Rulings R1-R4. No default changes; Phase 4 decides per client.
```
