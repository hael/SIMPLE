# Continuous view refinement on polar sections: review report (Phases 0-2)

Status, 2026-11-19: the run `polar_continuous_view` (arun, FWM-S116118, unattended) carried out Phases 0-2 of
`polar_continuous_view.md` on base `7e88a551e`, under rulings R1-R4. Phase 3 (the PorF arms and the verdicts) is
Hans's (R1); section 8 gives him the commands and the rule. The plan stays in `planned/`. Nothing is committed. The
plan's Progress section holds the detail; evidence paths are under `/Users/elmlundho/agent_work/arun/polar_continuous_view`.

## 1. What was learned

- **The B-spline sections are as good as the KB bank.** On the 4-blob phantom (shells 2-17, the 20 directions of the
  `nspace=20` spiral), the section's relative error is 5.1-8.3 % against the KB bank's 4.9-7.9 %, and the normal
  derivative's 5.6-7.7 % against the closed-form KB derivative's 9.4-12.7 %. T1 and T2 pass with margin
  (floors 2× KB, frozen in Phase 0). The gather is the 3-tap quadratic B-spline of `cartft_calc`, on references
  pre-weighted to reproduce the KB gather on the grid; a section costs what a KB section costs (153-160 µs against
  153-157 µs, 72 × 30 samples, one thread).
- **The five-variable cost is exact.** Its gradient matches central differences (worst 1.5e-4 against a tolerance ≥
  3e-3), its Gauss-Newton matrix is the Hessian at zero residual (≤ 2.6e-5 relative), and both checks fail on planted
  defects (a sign-flipped dV/dn, a flipped Gauss-Newton column).
- **The view is recovered, and LM is far cheaper here.** A noiseless particle 3° off the nearest of 1,000 cells,
  17.3° in-plane and (1.4, −0.6) px off is recovered to 0.055° (`lbfgsb`) and 0.060° (`lm`) in view (floor 0.25°);
  with the view off it keeps its 3.00°. Sections per solve on that fixture: 65 (`lbfgsb`) against 4 (`lm`). On the
  `cont_refine3D_1jxy` smoke runs (1,000 particles, 4 iterations, nspace = nspace_sub = 1,000): 53-56 sections per
  solve for `lbfgsb` against 15 (iteration 1) and 10 (iterations 2-4) for `lm`, 53 s against 35 s of wall time, and
  the same final resolution, FSC 0.143 / 0.5 = 3.17 / 3.60 Å, against 3.47 / 3.98 Å with the view off at the same
  nspace (30 s). A repeat gave the same counts and resolutions. These are smoke runs on a simulated set, not a verdict.
- **A base defect in the polar sigma2 (not changed; section 7).** `calc_sigma2_pftc` gets the stored-frame shift
  increment and applies it as the native shift. At a nonzero in-plane angle the residual is that of a wrongly placed
  shift: 1.17e-2 against the joint cost's 6.26e-3 on the recovery fixture. View-committed particles take the cost's
  own residual, so the H arms' sigma2 differs from A's and C's by this, not only by the view.

## 2. What changed, phase by phase

| Phase | Diff | What |
| --- | --- | --- |
| 0 | `review/phase_0.diff` | Measurement only (the plan): base builds and gates, the KB errors of T1/T2 by a scratch program (R2) and their frozen floors, the no-change case (63 identity outputs, exclusions), rulings R1-R3 in the plan |
| 1 | `review/phase_1.diff` | `cartft_calc%gather_ref`; `polarft_cont_calc` (new: sections, dV/dn, per-thread cache, the tilt chart, the cost); `polarft_calc%gen_sect_euclid_grad_at_angle` (the evaluator on supplied sections, its Gauss-Newton matrix) with `heap_vars%view_coeffs`; T1-T4 and their negative controls |
| 2 | `review/phase_2.diff` | The key `cont_view` (parameters, UI, commander refusals, nspace_sub = nspace); the 5-variable mode of `pftc_shsrch_grad`; the view branch and pose commit of `strategy3D_srch`; both representations, the commit, the final-section sigma2 and BENCH counts in the matcher; `polarft_cont_calc` in the builder; the recovery and refusal tests; the policy; this report |

`review/full_after_N.diff` holds the cumulative diff.

## 3. Decisions taken under delegation

All are recorded in Progress, under each phase.

- **Phase 0.** T1/T2 on the unmasked phantom (the masked one is not analytic) at shells 2-17 (analytic RMS ≥ 1e-6 of
  shell 2's), lp 2.7 Å; the KB reference of T2 is the closed-form derivative of the local KB interpolant;
  `uptime.log` joins the exclusions.
- **Phase 1.**
  1. Decision 1's "cubic B-spline" is `gather_packed_bspline3_window_grad`, the gather section 1 names (3 taps,
     quadratic).
  2. The cost's variables are the joint search's (shift x, shift y, rotation index) plus the tilts of
     T(w) = exp([wx, wy, 0]x) about the drawn reference's in-plane axes, so that grad(1:3) is T3's.
  3. The tilt columns are the exact derivatives of the section in that chart (left Jacobian of the exponential map).
  4. The Gauss-Newton matrix is a direct sum like `gen_joint_curv_at_angle`'s, not 15 angular-FFT columns.
  5. The new tester runs from the "continuous in-plane" suite: a suite of its own needs the test UI's list
     (`ui/simple_test/simple_test_ui_class.f90`), outside the table.
  6. Fixture and FD constants are copied into the new testers (no entity public for a tester; a shared fixture module
     would be outside the table).
- **Phase 2.**
  1. Under the key, `nspace` defaults to 1,000 (decision 8's "(1,000)") and the commander sets `nspace_sub = nspace`.
  2. The pose commit: `assign_ori` and its caller are outside the table, so the matcher calls
     `strategy3D_srch%commit_view` right after the search. It overwrites the Euler triplet with
     euler2m(assigned)·R_refᵀ·T(w)·R_ref and sets `proj` to the nearest grid direction. The convergence distances stay
     on the discrete cells (`assign_ori`'s own rule).
  3. The view pass reads the B-spline references into a calculator of its own (`read_cart_refvols_view`), not
     `build%cftc`. References there make every batch prepare Cartesian particle slots (`build_batch_particles3D`,
     outside the table); the first smoke runs stopped on exactly that.
  4. Sigma2: only a committed (improving) view solve takes the residual of its final evaluation; every other particle
     keeps the bank's, as before. `gen_sigma_contrib` was split into a loader and one shared body, so the bank path
     runs the same operations (the no-change case is byte-identical).
  5. A non-improving view solve stores the B-spline re-score of the drawn pose as `corr` (the contract's "B-spline
     scores only"); an invalid solve keeps the table's.
  6. `cont_view` with `trs=0` is refused (it would do nothing but force the grid).
  7. The recovery fixture takes the first tilt direction (30° steps) that keeps the drawn cell the truth's nearest
     direction (the plan's "off the nearest of 1,000 cells"); the 40° first try did not, and the fixture's own
     precondition caught it.

## 4. Deviations from the plan

- **Plain `refine3D` accepts the key** on `refine=prob_neigh` with one state and `objfun=euclid`. It is the stage
  program `refine3D_local` runs in process, and its workers run as `refine3D` with the key; only the router
  (`exec/simple_exec_refine3D.f90`, outside the table) can tell a user from the stage. `cont_prior`, another
  `refine3D_local` key, has the same standing. Every other program, route and setting is refused before workers
  launch (`scratch/logs/refusals.txt`: `refine=cont`, nspace_sub ≠ nspace, a bad value, two states,
  `refine3D_auto`), and other programs' parameter construction ignores the key (the master's reconstruction children
  inherit its command line).
- **The Gauss-Newton matrix by a direct sum** (Phase 1, decision 4 above), not the plan's "15 extra angular-FFT
  columns".
- **Rows of the table not used:** `simple_eul_prob_tab_neigh.f90` (with nspace_sub = nspace the subspace map is the
  identity; nothing to change), `commanders/test/simple_commanders_test_class.f90` (the suite).
- **The chart's sign.** dS/dwx = −k_y dV/dn in the code's chart (section 1 writes +k_y): a sign convention, verified
  by T4 and by the recovery.

## 5. Evidence against the baseline

- **Linux baseline-check portability.** The frozen KB errors were measured on macOS with GCC 16.
  Linux/GCC 15.2 reproduces the low shells but gives 0.079582 at shell 17 against the frozen
  0.079389 (0.24 % apart), failing the original 0.1 % reproduction check. That shell's analytic
  RMS is only 1.23e-6 of shell 2's: the original check demands absolute agreement below
  1e-10 of shell 2, despite the single-precision phantom and FFT. The reproduction check now
  adds one single-precision epsilon of shell 2's analytic RMS, converted to each shell's
  relative-error units. The allowance is capped at the frozen KB error, so the bank must
  still satisfy the frozen T1 ceiling. A zero-bank negative control must fail at every shell.
  The frozen T1/T2 accuracy limits, judged shells and production calculations are unchanged.
  Validation on Linux/GCC 15.2 Release: incremental build and installation passed; the full
  `unit_pftc_align2D3D` suite passed in reverse order; the fast gate passed all 15 entries in
  22.40 s, with 6,602/6,602 assertions passing in `unit_pftc_align2D3D`.
- **No-change case** (`cont_view=no`, Release, `nthr=1`): the 63 identity outputs byte-identical to Phase 0's base run
  after every phase (`scratch/logs/cmp_p1_prob_neigh_final.txt`, `cmp_p2_prob_neigh_final.txt`).
- **Gates:** fast gate 15/15 on Release and Debug; `fills` 48 rows, no differing entry, its table identical to the
  base's but for wall time; T4 with its negative controls.
- **Suites** (Release and Debug): `unit_pftc_align2D3D` 6471 assertions (base 6434), `unit_core` 1138,
  `unit_cart_align3D` 1178 (base 1176).
- **Smoke runs** (`scratch/runs/smoke_*`, repeated in `smoke2_*`): all three complete; the bank file holds exactly
  1,000 sections per half under the key.
- **T1-T4, the recovery, the sigma2 residual, the refusals:** section 1 and Progress.

## 6. Proposed skill text (not applied: skills are not edited in this run)

- **`simple-main-pftc`, Rules:** "`gen_sect_euclid_grad_at_angle` is `gen_raw_euclid_grad_at_angle` against a supplied
  section and its two view derivatives (`polarft_cont_calc`), with the gradient and the Gauss-Newton matrix in (shift
  x, shift y, rotation index, tilt x, tilt y); it reuses the thread's FFT buffers and `heap_vars%view_coeffs`.
  `gen_sigma_contrib_sect` is `gen_sigma_contrib` against a supplied section. `polarft_cont_calc` (`build%pcont`) holds
  pointers to `pftc` and a `cartft_calc`: it must be renewed with `pftc`, and its calculator must be built on
  `pftc`'s box. `set_view_ref` before a solve, on the solving thread." **Where to start:** add
  `simple_polarft_cont_calc.f90`. **Tests:** the polar continuous view tests run in the "continuous in-plane" suite.
- **`simple-main-cftc`, Rules:** "`gather_ref(state, iseven, q, val, grad)` is the point query of the polar view
  refinement: q on the native lattice, the value in a polar section's scale, the gradient per native unit."
- **`simple-refine3d`:** "`cont_view=lbfgsb|lm` (refine3D_local refine=prob_neigh, one state) refines the committed
  pose's view on B-spline sections after the draw: nspace_sub = nspace, the commit by `strategy3D_srch%commit_view`
  in the matcher after the search, the final section's sigma2 (`calc_sigma2_view`), the view pass's own
  `cartft_calc`. Policy: refine3D_policy.md, 'Continuous view refinement on polar sections'."

## 7. Left open for the maintainer

- **The PorF verdicts** (Phase 3, R1): section 8.
- **The base sigma2 frame (fixed 2026-11-19, after the run).** `calc_sigma2_pftc` (`sigma2/simple_euclid_sigma2.f90`)
  applied the matcher's stored-frame shift increment as the native shift of `gen_sigma_contrib`: on the recovery
  fixture (view off, in-plane 17.3°) a residual of 1.17e-2 against the joint cost's 6.26e-3. It now rotates the
  increment into the particle frame at the evaluated rotation. This changes every polar refinement's sigma2 (2D and
  3D) at a nonzero in-plane angle and shift; A, C and H all run on the fixed build.
- **`refine3D` user refusal:** kept as is (Hans, 2026-11-19).
- **Test suites: the existing ones (Hans, 2026-11-19).** T1-T4 run in `continuous in-plane`, the recovery and
  commit in `refine3D in-plane state`, the sigma2 frame in `sigma2 anchor`, the key in the parameters tester; no
  suite of its own. No registered high-level entry runs `cont_view`; add a case to one when the default switches.
- **The default switch** (decision 6): the key stays `no` until the PorF read; it becomes `no|yes` with the winning
  solver.

## 8. For Hans: the PorF arms and the verdict rule

Run each command twice per arm, on the build of this change, from `PorF_refine3D_local_tests` (the `local4_best`
command, `~/for_claude/PorF_refine3D_local_tests/local4_best/log_local`):

- **A** (baseline; nspace 20,000, nspace_sub 500 by default):
  `simple_exec prg=refine3D_local pgrp=c1 mskdiam=170 nparts=10 nthr=8 projfile=3_solve3D/tst.simple rec_backend=pcg refine=prob_neigh`
- **C** (control): A's command plus `nspace=1000 nspace_sub=1000`
- **H-lm:** A's command plus `nspace=1000 nspace_sub=1000 cont_view=lm`
- **H-lbfgsb:** A's command plus `nspace=1000 nspace_sub=1000 cont_view=lbfgsb`

(`nspace_sub` may be left out under `cont_view`: it is set equal to `nspace`, and a different value is refused.)

**What to compare.** For each run: the final FSC 0.143 resolution and the FSC at 4.22 Å; from
`BENCH/REFINE3D_BENCH_PART*.txt` at the last iteration, `match3D view sections per solve` (H only),
`match3D view evaluations per solve`, `match3D total thread-s` and `match3D orientation search`; peak RSS; wall time;
paired pose agreement with A (`porf_pair.py`); the maps.

**The verdict rule** (section 4 of the plan, pre-registered):

- s is A's two-run spread, in final FSC 0.143 and in FSC at 4.22 Å.
- **View refinement:** an H arm is a lead if its mean FSC 0.143 is within max(s, one shell) of A's and it beats C by
  more than s. If C is within s of A, the data cannot judge the view refinement: stop and report.
- **Solver:** the better FSC 0.143 wins when the two solvers differ by more than s; otherwise the solver with fewer
  sections per solve wins.
- With two runs per arm the result is a lead, not a finding.
- Section 7's sigma2 frame is fixed in the build all four arms use, so it no longer confounds H against A and C.
