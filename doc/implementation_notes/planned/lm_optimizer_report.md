# One Levenberg-Marquardt for the library: review report (Phases 0-4)

Status, 2026-11-18: the run `lm_optimizer` (arun, FWM-S116118, unattended) carried out Phases 0-4 of
`lm_optimizer.md` on base `aae467bc7`, under rulings R1-R10.

- **Phases 0-3** built the engine, the factory's `'lm'` and the clients' curvatures.
- **Phase 4** did the following:
  - timed both solvers quietly;
  - cached the polar joint curvature (R7);
  - fixed the factory LM's re-evaluation after callbacks (R10);
  - judged each client by R6, on its comparison fixture and on the production-shaped fixtures the tests already
    have.
- **The switches.** Two clients switch to `'lm'`: the polar joint (shift and in-plane) mode and `quality_logistic`.
  Three keep `'lbfgsb'`, because they are slower under `'lm'` on noisy data: the polar shift mode, `ftexp_shsrch` and
  the continuous CTF fit.
- **Still open.** The real-data reads are left to Cyril and Hans (R8, sections 7 and 8). The plan stays in
  `planned/` until they are in (R5).

Nothing is committed. The plan's Progress section holds the detail; evidence paths are under
`/Users/elmlundho/agent_work/arun/lm_optimizer`.

## 1. What was learned

- **One engine carries every LM client, the Cartesian route byte for byte** (`refine=cont`, Release and Debug,
  rechecked after every phase, Phase 4 included).
- **The step control as extracted did not suit every client** (the Phase 2 stop). The relative stops became options
  that take the factory's `factr`/`pgtol` (R3), and exact null coordinates take a zero step (R4).
- **Every client's Gauss-Newton curvature equals its Hessian at a zero residual** (2e-12 to 9e-5 relative). This
  includes the polar joint mode after R7's cache.
- **Noise-free fixtures flatter Gauss-Newton.**
  - On Phase 3's comparison fixtures (noise-free copies, zero residual), `'lm'` is a Newton method. It wins every
    client: 0.29-0.63 of L-BFGS-B's median time per solve, quiet, far beyond the A-A spread.
  - On noisy data the Gauss-Newton matrix lacks the residual's second-order term, and LM converges linearly to the
    same `factr` stop. On the production-shaped fixtures the tests already have, three clients were slower:
    - the polar shift mode in the `prob_neigh` case of `cont_refine3D_1jxy` and the fills bench;
    - `ftexp_shsrch` and the CTF fit in `simulated_workflow_6vxx`'s motion-correction and CTF steps.
  - The polar joint mode was no slower there. `quality_logistic` has the exact Hessian (Fisher plus ridge), so it
    is unaffected.
- **A defect in `opt_lm`, now fixed (R10).** After every accepted step, `opt_lm` re-evaluated whenever a callback
  existed. In the polar alternating search the callback changes the objective only when the rotation changes.
  - With the fix, the fill searches' evaluations per search under `'lm'` fall from 17.7 / 12.0 / 11.2 / 9.6 to
    11.5 / 6.5 / 6.1 / 5.3 (iterations 1-4; L-BFGS-B 7.2 / 6.5 / 6.6 / 5.5).
  - The shift mode's thread-second excess falls from 36 / 29 / 26 / 22 % to 34 / 20 / 17 / 13 %.
  - What remains is more evaluations in iteration 1, which fits Gauss-Newton on large residuals. Later it is more
    time per evaluation at equal counts, where the likely cause is the all-rotation scoring the callback does at
    each accepted step (not measured separately).
- **The CTF cost is several |sin| periods per 0.1 µm wide at 5 Å.** The continuous fit is a local refinement of the
  grid search.

## 2. What changed, phase by phase

| Phase | Diff | What |
| --- | --- | --- |
| 0 | `review/phase_0.diff` | Measurement: base builds and gates, identity cases and exclusions, per-client baselines, Cartesian outcome counts, the `cont_route` probe |
| 1 | `review/phase_1.diff` | `simple_lm_solver` (new), `spd_solve_rel`, `refine_joint`/`refine_shift` on the engine, `cont_route` retired, policy and two skills |
| 2 | `review/phase_2.diff` | `opt_lm` (new), `fdfcurvfun_8`, the factory's `'lm'`, L1-L6, R3 and R4 in the engine |
| 3 | `review/phase_3.diff` | Curvatures and solver names of the five factory clients with oracles and comparisons, the Guinier fit on the engine, the interim report |
| 4 | `review/phase_4.diff` | R7's joint-curvature cache (`load_ref_joint_curv_8`, its `polarft_core` allocation); the shift-mode curvature formed once per rotation; R10 (`opt_spec%callback_changed`, `opt_lm`, the polar callback, L4b, the polar re-evaluation count); median timing in the four client testers; `SHSRCH_SOLVER` split into the shift mode's (`'lbfgsb'`) and `SHSRCH_JOINT_SOLVER` (`'lm'`); `LOGISTIC_SOLVER = 'lm'`; the `simple-main-opt` skill; the plan (R5-R10, verdicts) and this report |

`review/full_after_N.diff` holds the cumulative diff. Phase 4's diff applies on top of Phases 0-3 (R5).

## 3. Decisions taken under delegation

Each is recorded, with its reason, in the plan's Progress entry of its phase.

- **Phases 0-3:**
  - the CTF baseline on `mini_stream`, the identity fixture, the exclusions and the probe method (Phase 0);
  - the `negligible_step` hook, the seed handling and N20's deletion, and export-built identity trees (Phase 1);
  - explicit `factr`/`pgtol` in L1-L6, and the L2 control at L2's tolerances (Phase 2);
  - each client's curvature, the polar log-form objective and the test seams (Phase 3).
- **Phase 4:**
  1. **R6's "the client's tester fixtures" includes the production-shaped fixtures the tests already have.**
     - The fixtures:
       - for the polar search, the `prob_neigh` case on `cont_refine3D_1jxy`'s fixture and the
         `polar_matching_bench` fills bench;
       - for `ftexp_shsrch` and the CTF fit, `simulated_workflow_6vxx`'s motion-correction and CTF steps;
       - for `quality_logistic`, its `cavg quality nano` suite.
     - The rule: a client switches when it is faster beyond the A-A spread on its comparison fixture and not slower
       beyond the spread on any of these.
     - The reason: the comparison fixtures are noise-free, and R6's intent is "good enough and fast" in use. Every
       negative is undone by one constant.
  2. **Each client has its own constant (decision 7).** The polar joint mode is its own client (its own `specify`),
     so the polar constant is split: `SHSRCH_SOLVER` for the shift mode and `SHSRCH_JOINT_SOLVER` for the joint mode.
  3. **R7's form.** Per solve, `load_ref_joint_curv_8` computes |REF|², Im(conj(REF)·REF′) and |REF′|², with one FFT
     batch pair for REF′. Each evaluation sums them with the analytic shift-phase derivative and the CTF² at the
     rotated position, linearly interpolated.
     - It equals Phase 3's FFT form at whole rotation indices up to aliasing (4th digit), and is positive
       semi-definite in between.
     - The shift-mode curvature, constant over the shift, is formed once per rotation of a solve (decision 6).
  4. **Medians.** The testers time repeated solves and report medians: pftc 20 × 3, CTF 5 × 8, `ftexp` 100,
     logistic 20.
     - "Faster than the A-A spread" means that the smaller of the two runs' gaps exceeds the larger of the two
       solvers' run-to-run differences.
     - On production-shaped fixtures the measure is thread-seconds of the affected stage, or wall time of a
       one-thread program run in a worker process.
  5. **R10's tester count.** It uses private profile counters of `pftc_shsrch_grad`, read through the existing
     test-only seam submodule (`shsrch_callback_profile`).
  6. **The cc shift-prior stop.** Under `'lm'` the shift mode stops on a cc objective with a shift prior. That
     combination is unreachable in production, and the shift mode stays on `'lbfgsb'`.

## 4. Deviations from the plan

- **Phase 2 stopped once** (R3, R4). Phase 1's `rg` exit holds as a word search, N20's stages i and j were deleted,
  and `dose series manifest` does not run the Guinier fit.
- **The R9 re-time of Phase 3's code printed means.** The verdicts use the median re-time taken after the tester
  change and R7.
- **Decision 8's real-data A/B in this run** was replaced by R5/R8.
- **R6 applied with production-shaped fixtures** (decision 1 above). By the comparison fixtures alone all five
  clients would have switched. A first B build did switch them all; it was measured, and three of the switches
  were reverted.
- **Timing under other people's load.** Some timings ran under someone else's load bursts (compiles, a load of 127
  at 12:46). The decisive pairs were repeated or agree run to run, and the uptimes are recorded beside each.
- **Files added in Phase 4:** `pftc/simple_polarft_core.f90` (R7), and Phase 4 on the living-descriptions row and
  the `opt` row (R10).
- **Stale wording outside the table, not edited:**
  - the `simple-main-pftc` skill: "the L-BFGS-B shift (or joint shift and in-plane) search" (the joint mode is now
    LM);
  - the `simple-cluster-cavgs-quality` skill: "its L-BFGS-B fit" (now LM).

## 5. Evidence against the baseline

- **Builds and gates.** Release and Debug without warnings. Fast gate on the final build: 15/15 on Release and Debug (`scratch/logs/fastgate_final_{release,debug}.log`).
- **Cartesian identity:** byte-identical on the declared outputs on the `'lbfgsb'` build after R10
  (`scratch/logs/cmp_r10_rel_cont.txt`) and on the final build (`scratch/logs/cmp_final_cont.txt`).
- **Polar identity, on the all-`'lbfgsb'` build after R7 and R10:** `prob_neigh` 63 byte-identical + 12 normalized
  (`scratch/logs/cmp_r10_rel_prob_neigh.txt`), and `fills` 48 rows with no differing entry. After the switch:
  48 rows, no differing entry (`scratch/runs/suites_final/fills/`).
- **`prob_neigh` on the final build** (result-changing, no gate): FSC 0.143 / 0.5 at 3.17 / 3.53 Å, cFAR 0.91, orientation overlap 0.998 and best-orientation distance 1.9388 deg as at the base; score 0.6396387 against 0.6396391, shift increment 0.013350 against 0.013340 px (`scratch/logs/cmp_final_prob_neigh.txt`, run `scratch/runs/r10_pn_A_1`).
- **Gates of the switched clients** on the final build: `unit_numerics` 2107 assertions, `unit_image` 982, `unit_pftc_align2D3D` 6337, `unit_cart_align3D` 1176 and `lib_cart_align3D` 52 pass (`scratch/runs/suites_final/`); `simulated_workflow_6vxx`: PASS, docking 0.9881, whole-volume 0.9067, FSC 4.23 Å (all-`'lbfgsb'` A: 0.9890, 0.9003, 4.30 Å);
  `mini_stream suite=6vxx`: 73 PASS, the base's six defocus errors exactly.

## 6. Verdicts (R6)

Comparison fixtures, quiet re-time (R9: one thread, two runs, load 2.5-2.9; `scratch/runs/p4_time_{1,2}/`), median
seconds per solve (timer resolution 1 µs). The polar shift-mode rows were taken before R10, which cut the
alternating rows' `'lm'` evaluations further (8.3 → 6.3 euclid, 9 → 7 cc).

| Client | Evaluations lbfgsb / lm | lm / lbfgsb time | Production-shaped fixture (A lbfgsb / B lm, two runs each) | Verdict |
| --- | --- | --- | --- | --- |
| polar shift mode (alternating, fixed; euclid, cc) | 10-44 / 4-9 | 0.38-0.60 | `prob_neigh` fill refinement thread-s, iteration 1: 12.72, 12.49 / 17.03, 16.77; iteration 4: 4.56, 4.46 / 5.07, 5.10 (with R10). Fills bench: with R10, one-thread wall time 1.53-4.20 times A's on every route, A-A under 1 % (`scratch/logs/r10_fills_ab.md`) | `'lbfgsb'` (slower in use) |
| polar joint mode (euclid, cc), with R7 | 19 / 4 | 0.29-0.34 | `prob_neigh` orientation search (holds it), before R10: 13.02, 11.49 / 11.44, 11.45 (it. 1) … 6.36, 6.16 / 5.12, 5.14 (it. 4): not slower | `'lm'` |
| `ftexp_shsrch` | 12.4 / 3.6 | 0.63 | `motion_correct` on the 16 simulated movies, one thread: 82.1, 96.9 / 122.7, 114.4 s | `'lbfgsb'` (slower in use) |
| `ctf_estimate_cost4Dcont` | 57.4 / 29.9 | 0.56 | `ctf_estimate` on the same: 42.1, 44.2 / 60.1, 68.9 s | `'lbfgsb'` (slower in use) |
| `quality_logistic` | 27 / 7 | 0.49 | `cavg quality nano` suite user time: 1.01, 1.08 / 1.07, 1.04 s (the fits are a small part) | `'lm'` |

No client has a failure under `'lm'` that L-BFGS-B does not have. Accuracy under `'lm'` meets every floor on all the
fixtures, and in `prob_neigh` the results agree with the base's to the printed digits.

## 7. For Cyril: motion correction and the CTF fit

- **What switched: nothing.**
  - Under `'lm'` both clients are faster on their noise-free tester fixtures. On `simulated_workflow_6vxx`'s noisy
    movies they took 18-50 % (motion correction) and 36-64 % (CTF estimation) more wall time.
  - Both stay on `'lbfgsb'`. Their `'lm'` paths are built and tested: oracles, comparisons, and `mini_stream` passing
    under `'lm'`.
- **The constant that switches each to `'lm'`:**
  - `FTEXP_SOLVER` in `src/main/image/simple_ftexp_shsrch.f90`;
  - `CTF_CONT_SOLVER` in `src/main/ctf/simple_ctf_estimate_cost.f90`.
- **What should not change beyond the fixtures' numbers under `'lm'`:**
  - shifts within 3.3e-4 px of the truth on the fixture (L-BFGS-B 1.4e-6; floor 0.05);
  - defocus within 0.002 µm and the angle within 0.04 deg;
  - on `mini_stream` the six defocus errors moved by at most 0.0001 µm.

  Expect real-data shifts, defocus and fit scores to agree to that order. The question for real data is time.
- **A/B commands.** Build A is the default (`'lbfgsb'`). Build B is the same snapshot with the constants switched
  to `'lm'`; this inverts R8's convention, because the default did not switch.

  ```
  cp -R SIMPLE SIMPLE_B   # the same snapshot; in SIMPLE_B:
  sed -i '' "s/FTEXP_SOLVER = 'lbfgsb'/FTEXP_SOLVER = 'lm'/" src/main/image/simple_ftexp_shsrch.f90
  sed -i '' "s/CTF_CONT_SOLVER = 'lbfgsb'/CTF_CONT_SOLVER = 'lm'/" src/main/ctf/simple_ctf_estimate_cost.f90
  # build each in its own tree (cmake -S <src> -B <tree> -DCMAKE_BUILD_TYPE=Release; cmake --build <tree> --target install);
  # then, with each build's environment, on the same movies, two runs per build:
  simple_exec prg=motion_correct projfile=<A|B>.simple nthr=<n>
  simple_exec prg=ctf_estimate   projfile=<A|B>.simple nthr=<n>
  ```

  Compare per-micrograph shifts, defocus, astigmatism and fit scores, and each program's wall time.

## 8. For Hans: the polar searches (PorF)

- **What switched.** The joint shift and in-plane mode of `pftc_shsrch_grad` now runs `'lm'`, with R7's cached
  curvature. Its users are the continuous in-plane refinement of refine3D/refine2D, `corrmat` and the 2D utilities.
  The shift mode (the tables' fill searches, the window and disk searches) stays on `'lbfgsb'`; its numbers are in
  section 6.
- **The constants:**
  - `SHSRCH_JOINT_SOLVER` in `src/main/pftc/simple_pftc_shsrch_grad.f90` switches the joint mode back
    (`'lm'` → `'lbfgsb'`);
  - `SHSRCH_SOLVER` switches the shift mode to `'lm'`, if you want to see it on PorF.
- **What should not change beyond the fixtures' numbers.** Resolution, orientation overlap and scores. On the
  `prob_neigh` case the joint mode changed none beyond the printed digits (FSC 0.143 / 0.5 at 3.17 / 3.53 Å, cFAR 0.91, orientation overlap 0.998 and best-orientation distance 1.9388 deg as at the base; score 0.6396387 against 0.6396391, shift increment 0.013350 against 0.013340 px (`scratch/logs/cmp_final_prob_neigh.txt`, run `scratch/runs/r10_pn_A_1`)). The cost of the continuous
  in-plane refinement should fall: its fixture solve takes about a third of L-BFGS-B's time.
- **A/B commands.**
  - Build B is the default. Build A is the same snapshot with
    `sed -i '' "s/SHSRCH_JOINT_SOLVER = 'lm'    /SHSRCH_JOINT_SOLVER = 'lbfgsb'/" src/main/pftc/simple_pftc_shsrch_grad.f90`,
    in its own tree.
  - Run the `local4_best` command twice per build, on the batch-capacity floor (`aae467bc7`).
  - Compare `match3D orientation search` and `match3D total thread-s` in `BENCH/REFINE3D_BENCH_PART1.txt` at
    iteration 8, the resolution and the orientation overlap.
  - For the shift mode as well, add a build C with `SHSRCH_SOLVER = 'lm'` and read `match3D fill refinement
    (thread-s)` and `match3D fill gradient evals per search`.

## 9. Left open for the maintainer

- **The real-data reads** (R8; sections 7 and 8).
- **A remedy for Gauss-Newton on large residuals,** if `'lm'` is wanted for the shift mode, `ftexp` or the CTF fit.
  Candidates are a curvature with the residual term, a looser `factr` under `'lm'`, or a cheaper discrete-angle
  callback. Each is a new decision.
- **Known limits of the engine:**
  - tiny but nonzero curvatures relative to the largest still fail the damped solve (R4 covers exact zeros);
  - the factory's coordinates are unscaled;
  - `'lm'` with a cc shift prior stops by design.
- **Stale skill wording outside the table** (section 4).

## 10. Commit message (for the maintainer, if the change is taken as one)

```
lm_optimizer: one Levenberg-Marquardt engine for the library

simple_lm_solver (lm_problem, lm_solve) owns the LM step control once: the Cartesian pose solve
(refine_joint, refine_shift) runs on it byte-identically, the factory gains 'lm' (opt_lm,
fdfcurvfun_8, L-BFGS-B's factr/pgtol as its stops, a re-evaluation only when a callback changed the
objective), the five factory clients gain verified Gauss-Newton curvatures (the polar joint one
cached per reference), and the dose-series Guinier fit runs on it. cont_route/shift_then_joint are
retired. By ruling R6 the polar joint search and quality_logistic default to 'lm'; the polar shift
search, ftexp_shsrch and the continuous CTF fit, slower under 'lm' on noisy data, keep 'lbfgsb'.
One private constant per client switches it. Rulings R1-R10.
```
