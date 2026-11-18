# Continuous view refinement on polar sections

Status, 2026-11-17: planned. Base: `e94681175`. Depends on: the batch-capacity floor in
`simple_matcher_memory.f90`, in the working tree (validated: no); LM in the optimizer factory,
Phases 1-2 of `lm_optimizer.md` (validated: no). If accepted, replaces sections 6.2-6.5 of
`refine3d_local_cartesian_priors_todo.md` (banner added then).

## 1. Why

PorF `refine3D_local` from `3_solve3D` (37,502 particles, C1, 8 iterations): `prob_neigh` gives
3.85/4.33 Å (FSC 0.143/0.5), `cont` 4.00/4.59 Å. 77 % of the paired poses that disagree by more
than 10° (5,300 of 6,894) include an in-plane error above 20°: the polar in-plane scan reaches
them, the Cartesian descent cannot. The polar price is memory, not compute: the 20,000-direction
bank takes 17.6 GiB per worker, while in iteration 8 refinement took 950 thread-s against 398 for
scoring. Aim: keep the polar search on a coarse bank of 1,000 directions, with nspace_sub = nspace
so that the bank holds only the coarse directions, and refine the drawn pose's view continuously.

**The view gradient.** A polar sample at (kx, ky) sits at loc = Rᵀ(kx, ky, 0). A tilt (ωx, ωy)
about the in-plane axes moves it along n = Rᵀẑ by ωx·ky − ωy·kx, so ∂S/∂ωx = ky·∂nV and
∂S/∂ωy = −kx·∂nV, where ∂nV = ∇V(loc)·n. This is one extra sampled quantity. It enters
`gen_raw_euclid_grad_at_angle` (`simple_polarft_corr.f90`) like the two shift columns, so the view
gradient comes out at every in-plane rotation from the same batched angular FFT.

**Checked at the base.**

- Every drawn 3D pose already gets a shift + in-plane solve (`refine_assignment_continuously`).
- The coarse cells of `prob_neigh` index into the bank (`subspace_inds`, `simple_builder.f90`), so
  the bank holds all nspace directions. With nspace_sub = nspace the bank holds only the coarse ones.
- `cart_refs_needed` already writes the Cartesian volumes beside the bank.
- Only the B-spline gather (`gather_packed_bspline3_window_grad`) has a gradient.
- Polar sigma2 is taken against the grid section `proj`.

## 2. Contract

- **Pipeline unchanged.** The fill, keyed draw, temperatures and weighted insertion stay as they
  are. The draw sees bank (KB) scores only. The refinement compares B-spline scores only: it
  commits on improvement over its own B-spline re-score of the drawn pose.
- **Ownership.** `cartft_calc` keeps the 3D B-spline references and gains one public point
  query, `gather_ref` (value and gradient). The new `polarft_cont_calc`
  (`pftc/simple_polarft_cont_calc.f90`) owns the polar side: the section and ∂nV at a rotation,
  the per-thread section cache, the tilt parametrisation and the solver's cost interface.
  `polarft_calc` keeps the particles and correlates supplied sections. `cartft_calc` and
  `polarft_calc` never import each other.
- **Committed pose.** Composed from matrices, Rz(ψ)·R_tilt·R_ref, never by pairing e1, e2 with
  a new e3. `proj` is the nearest grid direction. The sigma2 residual is the final evaluation's.

## 3. Decisions (2026-11-17, Hans)

1. **Interpolator.** Cubic B-spline on `cartft_calc`'s padded references for on-the-fly
   sections; the KB bank stays for coarse scoring. KB derivatives (the KB window's derivative is
   closed form) are the fallback if T1/T2 show the two scales disagree.
2. **Owner.** The solve of `refine_assignment_continuously` gains two view variables. No new
   strategy, route or program.
3. **Bound.** Tilts about the drawn reference's in-plane axes, bounded by
   |ωx|, |ωy| ≤ 2√(4π/nspace) rad, two grid spacings (12.8° at nspace 1,000; ruled "at least 2×").
4. **Objective.** Polar Euclidean only. The two view columns use the `+conjg` Friedel extension;
   their |S|² derivative columns are 2Re(S̄·∂S).
5. **Poses refined.** Only the committed pose; weighted-insertion extras keep the 3-variable
   solve. Refining several posterior samples is section 6, not this plan.
6. **Switch.** One key, `cont_view=no|lbfgsb|lm`, on `refine3D_local refine=prob_neigh` with one
   state only; refused elsewhere before workers launch. It stays `no` until the real-data read,
   and it becomes `no|yes` once decision 7's comparison has a winner.
7. **Optimizer: both, compared (ruled).** The joint object of `pftc_shsrch_grad` gains a
   5-variable mode that runs `'lbfgsb'` or `'lm'` from the factory (`lm_optimizer.md`).
   - **One cost interface.** `polarft_cont_calc` supplies f, the gradient and the Gauss-Newton
     matrix (15 extra angular-FFT columns); L-BFGS-B ignores the matrix.
   - **One chart.** (ωx, ωy, ωz, tx, ty) is an increment from the drawn pose by the exponential map,
     with ωz the in-plane angle; derivatives are taken in that chart.
   - **The case for LM, to be tested.** `cont`'s LM takes 3.3-9 passes per PorF solve, while the
     polar shift searches take 9-15 gradient evaluations. Here every evaluation that moves the view
     costs a section.
8. **Bank = coarse grid (ruled).** With `cont_view` on, nspace_sub is forced to nspace (1,000), and
   an explicit different value is refused. The coarse search scores the bank itself, and no
   reprojections are made beyond its 1,000 directions. The continuous view refinement takes the
   place of the fine neighbourhood step.

Defaults taken without a question: the key name.

**Not in this plan:**

- Multi-sample refinement (section 6).
- LM itself and its other clients (`lm_optimizer.md`).
- `refine=cont`, `simple_cartft_pose_opt*` and `simple_strategy3D_cont*`.
- solve3D and `refine3D_auto` stages.
- Multi-state, the cc objective, GPU offload; KB derivatives unless decision 1's fallback is
  ruled in.
- New programs, CMake targets or record slots.

## 4. Acceptance and stops

Unit tests use the 4-blob phantom of `simple_pftc_inplane_tester.f90`, whose transform is
analytic:

- **T1, section accuracy.** The B-spline section's per-shell relative RMS error against the
  analytic transform is ≤ 2× the KB bank's error at the same 20 rotations. The KB error is
  measured in Phase 0; the factor 2 is guessed.
- **T2, normal derivative.** ∂nV against the analytic ∂F/∂n, with the same rule.
- **T3, identity.** With the view frozen and the bank section supplied, f and grad(1:3) equal
  `gen_raw_euclid_grad_at_angle` within `IDENTITY_ATOL` (1e-4, an existing constant).
- **T4, gradient.** Central differences in ω (`FD_STEP`) over `PROBE_SHIFTS` × `PROBE_ROTOFFS`
  agree within max(`GRAD_FD_RTOL`·max|g|, `GRAD_FD_AFLOOR`) (existing constants). A sign-flipped
  ∂nV must fail. The Gauss-Newton matrix equals the central-difference Hessian at zero residual.

**PorF verdict, pre-registered.** Two runs per arm:

- **A:** nspace 20,000, nspace_sub 500, view off (the baseline, today's defaults).
- **C:** nspace = nspace_sub = 1,000, view off (the control).
- **H-lm** and **H-lbfgsb:** nspace = nspace_sub = 1,000, view on.

s is A's two-run spread in final FSC 0.143 and in FSC at 4.22 Å (measured, frozen after Phase 0).

- **View refinement.** An H arm is a lead if its mean FSC 0.143 is within max(s, one shell) of A
  and it beats C by more than s. If C is within s of A, the data cannot judge the view
  refinement: stop and report.
- **Solver.** The better FSC 0.143 wins when the two solvers differ by more than s. Otherwise the
  solver with fewer sections per solve wins.

With two runs per arm the result is a lead, not a finding.

Reported, not gated:

- alignment thread-seconds, peak RSS and wall time;
- paired pose agreement with A (`porf_pair.py`);
- the maps, which Hans reads.

**Stops.** T1 or T2 fails (the B-spline and bank scales disagree; ask whether to take decision 1's
KB fallback); sigma2 would need a parallel
path; a file outside the table; the UI refuses an arm command; PorF is unreachable (Phase 3
only).

**Gates no ruling waives.** The no-change case (with `cont_view=no`, byte-identical to its two
base runs at `nthr=1`), T4 with its negative control, the fast gate,
`polar_matching_bench suite=fills`.

## 5. Phases

| Phase | Kind | Content | Exit |
| --- | --- | --- | --- |
| 0 | measurement | Release + Debug builds with the floor and `lm_optimizer.md` Phases 1-2; fast gate, `unit_pftc_align2D3D`, `fills`; KB errors for T1/T2. No-change case: `refine3D_local refine=prob_neigh`, `nthr=1`, on the `cont_refine3D_1jxy` project, run twice with exclusions declared (stop if it cannot run). Arms A (the `local4_best` command) and C (the same plus `nspace=1000 nspace_sub=1000`) twice on the PorF workstation (Hans). | All recorded in Progress, s frozen |
| 1 | tests-only | `cartft_calc%gather_ref`; `polarft_cont_calc` (section + ∂nV at a rotation, cache, one cost interface for both solvers); the `polarft_calc` evaluator on supplied sections (gradient and Gauss-Newton matrix); per-thread buffers, no allocation in the loop | T1-T4 pass; new routines called only from testers; ns per section (B-spline vs KB) recorded |
| 2 | result-changing (`cont_view=lbfgsb\|lm`) | Key with nspace_sub forced to nspace, 5-variable joint mode, view branch, `polarft_cont_calc` in the builder, both representations read, `proj` update, final-section sigma2, sections-per-solve counter in BENCH, policy text | The no-change case and `fills` pass; with the view on, the bank holds exactly nspace reprojections and a different explicit nspace_sub is refused; a noiseless particle 3° in view, 17.3° in-plane and (1.4, −0.6) px off the nearest of 1,000 cells is recovered to ≤ 0.25° in view (guessed) by both solvers, and keeps 3° with the view off; sigma2 residual equals an independent recomputation within 1e-4; key refused off-route; sections per solve recorded for both solvers on that fixture; short smoke run per solver completes |
| 3 | result-changing | H-lm and H-lbfgsb twice each on PorF; verdicts; `polar_continuous_view_report.md` | Verdicts and report in Progress |

### File table

| File | What | Phases |
| --- | --- | --- |
| `cftc/simple_cartft_calc.f90`, `_tester` | `gather_ref` (public point query) | 1 |
| `pftc/simple_polarft_cont_calc.f90`, `simple_polarft_cont_calc_tester.f90` (new) | Sections, cache, cost interface; T1, T2, T4 | 1-2 |
| `pftc/simple_polarft_calc.f90`, `simple_polarft_corr.f90`, `simple_pftc_inplane_tester.f90` | Evaluator on supplied sections; T3 | 1 |
| `commanders/test/simple_commanders_test_class.f90` | Register the new tester in `unit_pftc_align2D3D` | 1 |
| `pftc/simple_pftc_shsrch_grad.f90` | 5-variable joint mode, solver by name | 2 |
| `strategies/search/simple_strategy3D_srch.f90`, `simple_strategy3D_inplane_tester.f90` | View branch, pose, `proj`, recovery | 2 |
| `simple_builder.f90` | Allocatable `polarft_cont_calc` beside `pftc` and `cftc` | 2 |
| `strategies/search/simple_strategy3D_matcher.f90`, `simple_matcher_refvol_utils.f90`, `sigma2/` (if needed) | Both representations, sigma2 | 2 |
| `params/simple_parameters{,_parse,_phases,_tester}.f90`, `ui/simple/simple_ui_refine3D.f90`, `commanders/simple/simple_commanders_refine3D.f90` | Key, validation, nspace_sub = nspace | 2 |
| `doc/policies/3D/refine3D_policy.md`, this note | Policy, Progress | 2-3 |

Paths under `src/main/` unless rooted.

## 6. Forward: refined posterior samples (not this plan)

The direction after this plan (Hans, 2026-11-17):

- **Sample.** Draw K samples from the tempered posterior over a bank of 500-1,000 directions,
  each with its own keyed draw. The draws are made, never an argmin.
- **Refine.** Refine each sample with this plan's 5-variable solve. The hook for this already
  exists: `refine_extra_pose_continuously` refines the weighted-insertion extras.
- **Merge.** Samples that converge within the insertion cone of one another
  (`insert_cone_from_fsc`) are one basin: they merge into one pose and their weights add.
- **Weight or draw.** For distinct basins, either each enters the reconstruction weighted by its
  refined score at the particle's own temperature, or the committed pose is a keyed draw among
  the basins.

The insertion stays within the drawn state, as the weighted insertion is now. The open question
is the weight-or-draw choice: the standing rule marginalises only inside the cone and leaves
hypotheses to the draw, so weighting distinct basins widens it and needs its own test. Refined
samples would replace the table's top-k extras and may let the bank shrink below 1,000. This
plan keeps the path open: the solve takes any drawn cell, and decision 5 only limits which cells
it is applied to.

## Progress

- 2026-11-17, session: plan written; premises checked at `e94681175`; nothing built.
