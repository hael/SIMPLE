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

Rulings (Hans, 2026-11-19, at launch; R1-R3 first applied in Phase 0):

- **R1. This run is Phases 0-2; PorF is Hans's.** Arms A and C of Phase 0 and all of Phase 3 run on
  the PorF workstation. s is not frozen in this run; nothing in Phases 0-2 depends on it. The review
  report of Phase 2 ends with the A, C, H-lm and H-lbfgsb commands and section 4's verdict rule for
  Hans. The plan stays in `planned/`.
- **R2. Where the KB errors of T1/T2 are measured.** Phase 0 edits no source file: the KB bank's
  per-shell errors come from a scratch program linked against the base build (or the first step of
  Phase 1's tester). The KB numbers and the T1/T2 floors (2× KB) are written into Progress and frozen
  before any B-spline section is evaluated against them.
- **R3. The file table.** Rewritten at launch with explicit paths so the fence reads it; the rows
  marked "(if needed)" (`simple_polarft_core.f90`, `simple_eul_prob_tab_neigh.f90`, the Phase 2
  wiring of `simple_polarft_calc.f90` / `simple_polarft_corr.f90`) are used only when the plan
  requires it, with the reason in Progress.
- **R4. T1 or T2 failing is a stop** (first applied in Phase 1, where both passed): stop and report the
  numbers so Hans can rule on decision 1's KB fallback; the run does not take the fallback itself.

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
  measured in Phase 0; the factor 2 is guessed. The setting (unmasked phantom, the nspace = 20 c1
  spiral, shells 2-17) and the per-shell floors of T1 and T2 are frozen in Progress, Phase 0 (R2).
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
| `cftc/simple_cartft_calc.f90` (and its `_tester`) | `gather_ref` (public point query) | 1 |
| `pftc/simple_polarft_cont_calc.f90` (new, and its `_tester`) | Sections, cache, cost interface; T1, T2, T4 | 1, 2 |
| `pftc/simple_polarft_calc.f90`, `simple_polarft_corr.f90`, `simple_pftc_inplane_tester.f90` | Evaluator on supplied sections; T3; Phase 2: wiring if needed | 1, 2 |
| `pftc/simple_polarft_core.f90` (if needed) | Per-thread buffers of the evaluator | 1, 2 |
| `commanders/test/simple_commanders_test_class.f90` | Register the new tester in `unit_pftc_align2D3D` | 1 |
| `pftc/simple_pftc_shsrch_grad.f90` | 5-variable joint mode, solver by name | 2 |
| `strategies/search/simple_strategy3D_srch.f90`, `simple_strategy3D_inplane_tester.f90` | View branch, pose, `proj`, recovery | 2 |
| `simple_builder.f90` | Allocatable `polarft_cont_calc` beside `pftc` and `cftc` | 2 |
| `strategies/search/simple_strategy3D_matcher.f90`, `simple_matcher_refvol_utils.f90` | Both representations, sigma2, sections-per-solve in BENCH | 2 |
| `strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` (if needed) | nspace_sub = nspace in the neighbourhood table | 2 |
| `sigma2/` (if needed) | Final-section sigma2 | 2 |
| `params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90`, `simple_parameters_tester.f90` | Key, validation, nspace_sub = nspace | 2 |
| `ui/simple/simple_ui_refine3D.f90`, `commanders/simple/simple_commanders_refine3D.f90` | Key in the UI, refusal before workers launch | 2 |
| `doc/policies/3D/refine3D_policy.md` | Policy | 2, 3 |
| `doc/implementation_notes/planned/polar_continuous_view.md` (this note), `polar_continuous_view_report.md` | Progress; the report | 0, 1, 2, 3 |

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
- **2026-11-19, Phase 0 (measurement), arun `polar_continuous_view`.** FWM-S116118: macOS 26.7.3
  (25G342), Apple M2 Ultra, 24 cores, 192 GiB; Homebrew GCC 16.2.0 (`gfortran-16`), cmake 4.2.1;
  builds at `-j12`. Base `7e88a551e` (HEAD of the copy; the floor `aae467bc7` and LM Phases 0-4 are in
  it), exported with `git archive` to `scratch/keep/base_src`. Load 5.0-10.1 during the session, from
  OS and IT daemons (backup client, media analysis, endpoint security), no one else's compute jobs;
  free disk 5.9 TiB. Paths below are under `/Users/elmlundho/agent_work/arun/polar_continuous_view`.
  Files changed: this note only (rulings R1-R3 recorded in section 3, T1's floor origin in section 4).
  - **Builds and gate.** Release and Debug base trees `scratch/keep/base_release`, `scratch/keep/base_debug`
    (`scratch/logs/build_base_{release,debug}.log`; no compiler warning, only the fortls registration
    notice of `CLAUDE_CLI=/usr/bin/false`). Fast gate 15/15 on both: Release 16.7 s, Debug 51.8 s
    (`scratch/logs/fastgate_base_{release,debug}.log`).
  - **Suites.** `unit_pftc_align2D3D` (1 thread, `OMP_THREAD_LIMIT=4`): 15/15 groups, 6434 assertions,
    on Release (11.5 s) and Debug (31.8 s) (`scratch/runs/suites_base_{release,debug}/`).
    `polar_matching_bench suite=fills nthr=8` (Release): 48 rows, no differing entry in any column,
    331.5 s, load 8.7-8.8 (`scratch/runs/suites_base_release/fills.log`, `scratch/logs/fills_base_release.tsv`).
  - **KB errors for T1/T2 (R2: a scratch program, no source changed).** `scratch/keep/kb_err/kb_bank_err.f90`
    (built against the Release base by `scratch/keep/kb_err/build.sh`; logs `scratch/keep/kb_err/run_*/out.log`).
    Setting, which Phase 1's tester reproduces: the 4-blob phantom of `simple_pftc_inplane_tester.f90`
    (box 64, 1.3 Å, its `CTRS`, `SIGMAS`, `AMPS`), **unmasked**; `lp` 2.7 Å (kfromto 2-31, pftsz 72); the
    20 rotations are the `nspace=20` c1 reference spiral (`build_refspiral`, e3 = 0; Euler triples in the
    log); the bank is `fproject_polar_batch_opt` on the 2× padded, expanded transform with `mirr_proj`
    off, exactly as `vol_pad2ref_pfts_opt` makes it (a re-implementation of its gather reproduces every
    bank value exactly). Analytic transform at q = Rᵀ(h, k_c, 0) (box lattice, `get_polar_coords`):
    F(q) = N⁻³ Σ_b A_b (2πσ_b²)^{3/2} exp(−2π²σ_b²|q|²/N²) exp(−2πi q·c_b/N), phase origin at voxel N/2+1;
    ∂F/∂n = Σ_b (−4π²σ_b²(q·n)/N² − 2πi c_b·n/N) F_b, per box voxel, n = Rᵀẑ. The KB normal derivative
    is the closed-form derivative of the normalized local KB interpolant (window at the point's own
    `nint`, `apod_fast_device`'s polynomial, double precision; decision 1's fallback); central differences
    of the single-precision gather agree at steps 0.03-0.1 voxel and are biased at smaller steps by the
    kernel's clamp at |x| = 1.5. Error per shell: sqrt(Σ|S − F|² / Σ|F|²) over 20 × 72 samples.
    - **Judged shells: 2-17**, where the phantom's analytic shell RMS is ≥ 1e-6 of shell 2's. Beyond,
      the transform falls to the single-precision rounding of the low shells (shell 18: 2.7e-7; the KB
      error rises to 0.11 at 19, 0.31 at 20, 1.9 at 21) and a relative error says nothing about
      interpolation. The fixture's 8 Å band (shells 2-10) gives the same numbers.
    - The soft-masked phantom (the fixture's mask) is not used: mask truncation dominates above
      shell 10 (relative error 0.13 at 11, 1.0 at 14; `run_lp2.7_m1_d0`), and the transform is analytic only
      without it.
    - **Frozen floors (origin: this measurement × the plan's factor 2):** T1 passes when the B-spline
      section's error is ≤ the T1 floor at every judged shell, T2 when ∂nV's is ≤ the T2 floor.

      | Shell | Å | analytic RMS / shell 2 | KB section error | **T1 floor** | KB ∂nV error | **T2 floor** |
      | --- | --- | --- | --- | --- | --- | --- |
      | 2 | 41.60 | 1.00e+00 | 0.0503 | 0.1007 | 0.1268 | 0.2536 |
      | 3 | 27.73 | 7.36e-01 | 0.0673 | 0.1347 | 0.0988 | 0.1977 |
      | 4 | 20.80 | 5.87e-01 | 0.0546 | 0.1092 | 0.0936 | 0.1872 |
      | 5 | 16.64 | 3.19e-01 | 0.0490 | 0.0980 | 0.0985 | 0.1971 |
      | 6 | 13.87 | 1.51e-01 | 0.0521 | 0.1041 | 0.0978 | 0.1955 |
      | 7 | 11.89 | 6.87e-02 | 0.0550 | 0.1100 | 0.1028 | 0.2056 |
      | 8 | 10.40 | 3.09e-02 | 0.0578 | 0.1156 | 0.1063 | 0.2126 |
      | 9 | 9.24 | 1.37e-02 | 0.0559 | 0.1118 | 0.0977 | 0.1955 |
      | 10 | 8.32 | 5.47e-03 | 0.0545 | 0.1089 | 0.1009 | 0.2018 |
      | 11 | 7.56 | 1.96e-03 | 0.0600 | 0.1200 | 0.1125 | 0.2251 |
      | 12 | 6.93 | 6.81e-04 | 0.0646 | 0.1292 | 0.1075 | 0.2151 |
      | 13 | 6.40 | 2.31e-04 | 0.0626 | 0.1252 | 0.1055 | 0.2110 |
      | 14 | 5.94 | 7.12e-05 | 0.0638 | 0.1275 | 0.1002 | 0.2004 |
      | 15 | 5.55 | 1.99e-05 | 0.0691 | 0.1382 | 0.1073 | 0.2146 |
      | 16 | 5.20 | 5.14e-06 | 0.0737 | 0.1474 | 0.1171 | 0.2341 |
      | 17 | 4.89 | 1.23e-06 | 0.0794 | 0.1588 | 0.1158 | 0.2316 |

      The KB section error is nearly flat over shells (5-8 %): the polar projector applies no gridding
      correction, so the kernel's real-space roll-off scales each blob by its distance from the centre.
  - **No-change case.** The identity fixture was captured again at this base by
    `scratch/capture_fixture.sh` (`lm_optimizer`'s recipe: `simple_test_exec test=cont_refine3D_1jxy
    nthr=1` on the Release base, stopped at its first iteration; version 1 of stage a's project is
    `enter_stage`'s copy of the seeded project). Fixture: `scratch/keep/identity_fixture/seeded_cont_gate.simple`
    (md5 `398eac55…`; byte for byte `lm_optimizer`'s apart from the stack's absolute path) with the maps
    and stack of `test_cont_refine3D_1jxy_96241/`. The case, `scratch/run_identity.sh <tree> <dir>
    prob_neigh`: `refine3D_local refine=prob_neigh vol1=1JYX.mrc mkdir=no pgrp=c1 mskdiam=120 nsample=500
    maxits=4 nthr=1 nparts=1`, `SIMPLE_SEED=20260923`, a fresh directory. It runs at `nthr=1`. Two runs,
    280.6 s each, peak RSS 4.26 GB (nspace 20,000: each half of the reprojection model is 1.34 GB);
    final FSC 0.143 / 0.5 3.23 / 3.60 Å. `scratch/compare_runs.sh` (`scratch/logs/cmp_base_prob_neigh.txt`;
    md5 lists `scratch/logs/md5_id_base_prob_neigh_{1,2}.txt`): 63 files byte-identical, the same set
    as `lm_optimizer`'s Phase 0.
    - **Identity outputs:** those 63 files: every `.mrc` (start, part, even/odd, unfiltered, NU-filtered
      and NU local-resolution, trailing, `rec_final*`, masks), `fsc*.bin`, the FSC `.eps`,
      `RESOLUTION_*`, `simple_stats.txt`, `sigma2_state.bin`, `trailrec_state01.txt`,
      `reprojection_model_{even,odd}.bin`, the `*_pcg_support.txt`, the two `.jpg`, the `*_FINISHED`
      markers, `EXIT_CODE_JOB_1`, `startvol_state01.mrc.local_import`, and the text dumps of the
      `ptcl3D` and `stk` segments (`print_project_field`).
    - **Declared exclusions** (fixed here, never extended): `cont_gate.simple` itself (absolute paths;
      judged by its `ptcl3D`/`stk` dumps, bytes, and its `out` and `projinfo` dumps, normalized);
      `run.log`, `BENCH/*`, `*.provenance` and `metadata.json`, judged after normalization (the run
      directory replaced, the TIMINGS, SPLIT TIMINGS and THREAD-SECONDS sections of BENCH dropped, lines
      with a wall time, RSS, PID or byte count dropped; all 12 NORM_SAME here); `EXIT_CODE_JOB_1.job`,
      `time.log` and `uptime.log` (host, PID, timing only; `uptime.log` matched here only because the two
      runs started together).
    - The first run is the "before" of every later identity check:
      `scratch/keep/identity_base/id_base_prob_neigh_1`.
  - **PorF.** Hans: A and C on the PorF workstation; s not yet frozen (R1).
  - **Exit.** Builds, gate, `unit_pftc_align2D3D`, `fills`, the KB errors and the no-change case are
    recorded above; s is Hans's (R1) and nothing in Phases 1-2 depends on it.
  - **Decisions under delegation.**
    1. T1/T2 are judged on the unmasked phantom (the masked one is not analytic) at shells 2-17 (the
       analytic RMS ≥ 1e-6 of shell 2's), at `lp` 2.7 Å over the 20 directions of the `nspace=20` spiral.
    2. The KB reference for T2 is the closed-form derivative of the local KB interpolant, not a finite
       difference of the bank (step-dependent at the kernel's clamp).
    3. `uptime.log` (written by this run's script) joins the exclusions as host timing; no other file was
       added.
    4. `unit_pftc_align2D3D` ran on both trees, `fills` on Release only (as in `lm_optimizer`).
  - What it means, in short: the base builds and passes its gates, the no-change case is deterministic
    at one thread with 63 identity outputs, and the B-spline section of Phase 1 must stay within twice
    the KB bank's 5-8 % section error and 9-13 % normal-derivative error at shells 2-17.
- **2026-11-19, Phase 1 (tests-only), arun `polar_continuous_view`.** FWM-S116118, toolchain as in Phase 0; base
  `7e88a551e`. Paths below are under `/Users/elmlundho/agent_work/arun/polar_continuous_view`. Trees `scratch/p1_release`
  and `scratch/p1_debug`, built from an export of the working tree without the repository metadata
  (`scratch/build_export.sh`, as the base was exported); no compiler warning in either
  (`scratch/logs/build_p1_release_7.log`, `build_p1_debug_3.log`).
  Files changed: `cftc/simple_cartft_calc.f90` (`gather_ref`, the public point query: B-spline value and gradient per
  native lattice unit), `cftc/simple_cartft_calc_tester.f90` (`test_gather_ref`), `pftc/simple_polarft_calc.f90` (the
  interface of `gen_sect_euclid_grad_at_angle`, `heap_vars%view_coeffs`), `pftc/simple_polarft_core.f90` (that buffer's
  allocation and release), `pftc/simple_polarft_corr.f90` (`gen_sect_euclid_grad_at_angle`, the evaluator on supplied
  sections), `pftc/simple_polarft_cont_calc.f90` (new), `pftc/simple_polarft_cont_calc_tester.f90` (new: T1, T2, T4,
  threads, cost), `pftc/simple_pftc_inplane_tester.f90` (T3; runs the new tester), this note. All are in Phase 1's rows.
  `simple_polarft_core.f90` is an "(if needed)" row (R3), used because the evaluator's two extra angular-coefficient
  columns are per-thread scratch, which the folder's rule puts in `heap_vars`, allocated in `new`. The row of
  `commanders/test/simple_commanders_test_class.f90` was not used (decision 6).
  - **What was built.** `polarft_cont_calc` holds pointers to `pftc` and `cftc`, the polar sample positions and one
    section cache per thread of the team at `new`. `set_view_ref(state, iseven, R_ref)` sets the calling thread's
    drawn reference; `eval(iptcl, x, f, grad, curv)` is the one cost interface of both solvers at x = (shift x,
    shift y, rotation index, tilt x, tilt y), the Gauss-Newton matrix optional (L-BFGS-B leaves it out); a section is
    cut only when the tilts change (`get_nsections` counts them). `section(state, iseven, F, S, dV/dn)` is the point
    query of T1/T2; `view_frame` gives F = exp([wx, wy, 0]x) R_ref. The evaluator reuses the thread's FFT buffers and
    plans (`cmat_joint_many`/`plan_fwd3_many`, `cmat2_many`/`plan_fwd1_many`, `crmat1_many`/`plan_mem_r2c_many`),
    with no allocation per call: per evaluation one 3-column batch (as `gen_raw_euclid_grad_at_angle`), two
    one-column batches and three real batches of nk transforms; the matrix adds two batches (the ring derivative) and
    an O(pftsz·nk) sum.
  - **Exit 1, T1-T4 pass.** `unit_pftc_align2D3D`, Release and Debug: 15/15 groups, 6457 assertions (Phase 0: 6434;
    T3 adds 3, the new tester 20), `scratch/runs/p1_final/pftc_{release,debug}.log`.
    - **T1 and T2** (`test_section_accuracy`, the setting frozen in Phase 0; the recomputed KB errors reproduce
      Phase 0's within 1e-3 relative, an assertion). Every judged shell passes:

      | Shell | KB section | B-spline section | T1 floor | KB dV/dn | B-spline dV/dn | T2 floor |
      | --- | --- | --- | --- | --- | --- | --- |
      | 2 | 0.0503 | 0.0535 | 0.1007 | 0.1268 | 0.0774 | 0.2536 |
      | 3 | 0.0673 | 0.0706 | 0.1347 | 0.0988 | 0.0724 | 0.1977 |
      | 4 | 0.0546 | 0.0573 | 0.1092 | 0.0936 | 0.0671 | 0.1872 |
      | 5 | 0.0490 | 0.0513 | 0.0980 | 0.0985 | 0.0629 | 0.1971 |
      | 6 | 0.0521 | 0.0545 | 0.1041 | 0.0978 | 0.0585 | 0.1955 |
      | 7 | 0.0550 | 0.0576 | 0.1100 | 0.1028 | 0.0599 | 0.2056 |
      | 8 | 0.0578 | 0.0608 | 0.1156 | 0.1063 | 0.0637 | 0.2126 |
      | 9 | 0.0559 | 0.0586 | 0.1118 | 0.0977 | 0.0610 | 0.1955 |
      | 10 | 0.0545 | 0.0570 | 0.1089 | 0.1009 | 0.0564 | 0.2018 |
      | 11 | 0.0600 | 0.0630 | 0.1200 | 0.1125 | 0.0610 | 0.2251 |
      | 12 | 0.0646 | 0.0677 | 0.1292 | 0.1075 | 0.0628 | 0.2151 |
      | 13 | 0.0626 | 0.0654 | 0.1252 | 0.1055 | 0.0590 | 0.2110 |
      | 14 | 0.0638 | 0.0664 | 0.1275 | 0.1002 | 0.0592 | 0.2004 |
      | 15 | 0.0691 | 0.0724 | 0.1382 | 0.1073 | 0.0628 | 0.2146 |
      | 16 | 0.0737 | 0.0779 | 0.1474 | 0.1171 | 0.0650 | 0.2341 |
      | 17 | 0.0794 | 0.0830 | 0.1588 | 0.1158 | 0.0705 | 0.2316 |
    - **T3** (`test_sect_evaluator_identity` in `simple_pftc_inplane_tester.f90`): with the bank reference supplied and
      zero view columns, f and grad(1:3) differ from `gen_raw_euclid_grad_at_angle` by at most 2.4e-8 and 1.0e-8
      (Release; 0 in Debug) and the Gauss-Newton block from `gen_joint_curv_at_angle` by 1.2e-16 relative, over the 3
      `PROBE_SHIFTS` at every 7th grid angle and the 4 `PROBE_ROTOFFS` beside each (`IDENTITY_ATOL` 1e-4).
    - **T4** (`test_view_gradient_and_curvature`, two cases: an even particle with unit sigma2, and an odd one whose
      half's sigma2 rises over the shells): the particle is the B-spline model at a truth pose (tilts (0.06, −0.04)
      rad, rotation index 11, shift (1.4, −0.6) px; loss there −1.1e-7 and 1.2e-8). The gradient matches central
      differences (`FD_STEP` 1e-3, all five variables) at the 12 `PROBE_SHIFTS` × `PROBE_ROTOFFS` poses with tilts
      offset by (0.025, −0.03) rad: worst error 9.9e-5 and 1.5e-4 (Release; 1.1e-4 and 1.4e-4 Debug), against
      max(`GRAD_FD_RTOL`·max|g|, `GRAD_FD_AFLOOR`) ≥ 3e-3. The Gauss-Newton matrix at the truth equals the
      central-difference Hessian of the analytic gradient within 1.4e-5 and 2.6e-5 relative (rule of `test_lm_curvature_oracles`: ≤ `GRAD_FD_RTOL`), also on the view
      rows alone, and is symmetric positive semi-definite.
    - **Negative controls.** In the tester: the gradient with dV/dn's sign flipped (grad(4:5) negated) fails the FD
      check (its error is 200× the tolerance: the largest gradient component is a tilt), and the matrix with the
      view-shift and view-rotation coupling flipped fails the Hessian check (4.2e-2 > 1e-2), in both cases. Planted in the code (in the scratch export only, then removed by
      re-exporting; `scratch/runs/mut{1,2}/`): M1, the sign of dV/dn in the tilt columns of `cut`, fails the gradient
      check (1.34) and the Hessian check (2.0); M2, the sign of one view column of J in the Gauss-Newton sum, fails only
      the Hessian check (0.37). After removal the suite passes again (`scratch/runs/p1_suite/cont2.log`).
    - **Threads** (`test_threads`): a team of four, each thread with its own drawn reference, gives f, the gradient
      and the matrix bit for bit as one thread; each thread's probes share one tilt, so the team cuts 4 sections where
      the serial pass, whose tilt changes at every probe, cuts 12 (the cache's reuse is asserted with a full team).
  - **Exit 2, new routines called only from testers.** `rg` outside `*_tester.f90` finds `gather_ref`,
    `gen_sect_euclid_grad_at_angle`, `view_coeffs` and `polarft_cont_calc` only at their definitions and in
    `polarft_cont_calc` itself, which only its tester and the in-plane tester use.
  - **Exit 3, ns per section** (one thread, Release, `test_section_accuracy`'s 200 sections each way, three runs, no
    job of mine): B-spline section with its normal derivative 153-160 µs, KB bank section (`vol_pad2ref_pfts_opt`
    over the 20 directions) 153-157 µs, both 72 × 30 samples (about 71-74 ns per sample), at load 2.8-3.1 from OS
    daemons, on the code before the review's guards (`scratch/runs/p1_suite`, recorded in the journal); repeated on the
    final code at load 4.2-4.5 (OS daemons): 153-154 µs both (`scratch/runs/p1_timing/timing.log`).
    The two cost the same; the cut of `eval` (two tilt columns instead of dV/dn) has the same gather.
  - **Gates.** No-change case on the final Release tree against Phase 0's base run: the 63 identity outputs
    byte-identical, every excluded file NORM_SAME or in the declared exclusions
    (`scratch/logs/cmp_p1_prob_neigh_final.txt`; 268 s; also after the first version of the code,
    `scratch/logs/cmp_p1_prob_neigh.txt`). `fills` (Release, 8 threads): 48 rows, no differing entry,
    the table identical to the base's but for wall time (`scratch/logs/fills_p1_release.tsv`). Fast gate 15/15:
    Release 17.7 s, Debug 53.2 s (`scratch/logs/fastgate_p1_{release,debug}.log`; its one warning names the two new,
    untracked sources). `unit_cart_align3D`: 1178 assertions, Release and Debug (1176 + `test_gather_ref`'s 2; the
    error lines it prints are its forked negative tests, printed by the base too).
  - **Decisions under delegation.**
    1. Decision 1's "cubic B-spline" is `cartft_calc`'s gather `gather_packed_bspline3_window_grad`, the one section 1
       names: a 3-tap quadratic B-spline (C¹; "3" counts the taps) on references pre-weighted so that on the grid it
       reproduces the KB gather. Used as is; T1 and T2 pass on it.
    2. The chart (decision 7): the variable order is the joint shift search's, (shift x, shift y, rotation index, tilt x,
       tilt y), so that grad(1:3) is T3's and the in-plane variable stays in rotation-index units; the tilts are the
       components of w in T(w) = exp([wx, wy, 0]x), radians, about the drawn reference's in-plane axes (decision 3);
       the pose is Rz(ψ) T(w) R_ref (section 2).
    3. The tilt columns are the exact derivatives of the section in that chart: dS/dw_a = −g·(v_a × q), g = F∇V the
       B-spline gradient in the section frame, v_a = J_l(w) e_a the left Jacobian of the exponential map; at w = 0 they
       are −k_y dV/dn and k_x dV/dn (section 1's derivation, the sign of the chart). The in-plane part of v_a (second
       order in w) comes from the same gather, so no ring derivative couples the tilts to the rotation index.
    4. The Gauss-Newton matrix is the direct sum of `gen_joint_curv_at_angle` extended by the two view columns (CTF²
       interpolated at the rotated position as there), not 15 extra angular-FFT columns: the solver needs it at one
       rotation, and the sum costs O(pftsz·nk) after one ring derivative.
    5. The section's power along the ring is transformed in each evaluation by `memoize_one_ref`'s r2c route, so a
       supplied bank section reproduces `ft_ref2`; the view columns' power terms are 2 Re(conj(S) dS) by the same route.
    6. The new tester runs from `run_all_pftc_inplane_tests` (the "continuous in-plane" suite of
       `unit_pftc_align2D3D`), not as a suite of its own: a new suite name must also be listed in the test UI's
       `suite=` help (`ui/simple_test/simple_test_ui_class.f90`, `check_test_registry`), a file outside the table. The
       precedent is `simple_ctf_tester`, which runs `run_all_ctf_estimate_cost_tests`. The report proposes the
       dedicated suite.
    7. The tester carries its own copy of the phantom and of the in-plane tester's FD constants (`FD_STEP`,
       `GRAD_FD_RTOL`, `GRAD_FD_AFLOOR`, `PROBE_*`), with the source named: the testing rules make no entity public
       for a tester, and a shared fixture module would be in `commanders/test/`, outside the table.
    8. T4's FD check covers all five variables, and the Hessian rule is also applied to the view rows alone (stricter
       than the plan's).
    9. After an independent review of the diff (no correctness finding), guards were added: `eval` stops when the
       polar grid changed since `new` or when the drawn reference's half is not the particle's, and `new` keeps
       pointers only to `intent(inout), target` calculators. Not tested: a CTF (`set_ptcl_ctfmat` is a private test
       hook of pftc's submodule testers); the CTF terms are those of `gen_raw_euclid_grad_at_angle` and
       `gen_joint_curv_at_angle`, extended to the view columns by the same expressions.
  - **For Phase 2.** `cartft_calc` must be made on the box of `pftc`'s polar grid (`gather_ref` returns the scale of
    its own box; `polarft_calc` exposes no box to check), and `polarft_cont_calc` renewed whenever `pftc` is. The
    committed pose is Rz(ψ)·`view_frame`(R_ref, w), with this phase's chart sign (dS/dwx = −k_y dV/dn).
  - **Finding (tooling).** `rsync -a` of the working tree keeps the files' times, so restoring a file in an export
    left its object file newer than the restored source and make did not rebuild it (seen after M2's removal).
    `scratch/build_export.sh` now copies by checksum without times.
  - What it means, in short: the B-spline sections cut from `cartft_calc`'s references are as accurate as the KB bank
    (sections 5-8 %, KB 5-8 %; normal derivative 6-8 %, KB's 9-13 %) and cost the same per section; the five-variable
    cost has an exact gradient and a Gauss-Newton matrix that is the Hessian at zero residual; nothing in production
    calls it, and the no-change case is byte-identical.
- **2026-11-19, Phase 2 (result-changing under `cont_view=lbfgsb|lm`), arun `polar_continuous_view`.** FWM-S116118,
  toolchain as in Phase 0; base `7e88a551e`. Paths below are under `/Users/elmlundho/agent_work/arun/polar_continuous_view`.
  Trees `scratch/p1_release`, `scratch/p1_debug` (rebuilt from the exported working tree, `scratch/build_export.sh`; no
  compiler warning, `scratch/logs/build_p2_{release_5,debug_1}.log` and later). R4 recorded in section 3.
  Files changed (all in Phase 2's rows): `params/simple_parameters.f90` (`cont_view`, `l_cont_view`),
  `simple_parameters_parse.f90` (binding), `simple_parameters_phases.f90` (value, route, nspace_sub = nspace),
  `simple_parameters_tester.f90` (`test_cont_view`); `ui/simple/simple_ui_refine3D.f90` (the input);
  `commanders/simple/simple_commanders_refine3D.f90` (nspace 1,000 and nspace_sub = nspace under the key, refusals);
  `simple_builder.f90` (`pcont`); `pftc/simple_pftc_shsrch_grad.f90` (`SHSRCH_JOINT_VIEW`, `new_joint_view`,
  `joint_view_terms`, `minimize_joint` for any dimension, `xsol`); `pftc/simple_polarft_calc.f90`, `_corr.f90`
  (`gen_sigma_contrib_sect` beside `gen_sigma_contrib`, one body `sigma_contrib_of_loaded_ref`: the wiring row);
  `pftc/simple_polarft_cont_calc.f90` (counters, `sigma_contrib`, `ensure_section`);
  `strategies/search/simple_strategy3D_srch.f90` (`view_optimizer`, `refine_assignment_view`, `commit_view`,
  `get_view_solution`), `simple_strategy3D_inplane_tester.f90` (`test_view_recovery`), `simple_strategy3D_matcher.f90`
  (both representations, the commit after the search, the final-section sigma2, BENCH counts, teardown),
  `simple_matcher_refvol_utils.f90` (`read_cart_refvols_view`, `cart_refs_needed`); `sigma2/simple_euclid_sigma2.f90`
  (`calc_sigma2_view`: the "(if needed)" row, used because the record is written there); `doc/policies/3D/refine3D_policy.md`;
  this note and `polar_continuous_view_report.md`. Not used: `simple_eul_prob_tab_neigh.f90` (nspace_sub = nspace
  needs no change there: the subspace map becomes the identity) and `simple_polarft_core.f90`.
  - **What was built.** The key `cont_view=no|lbfgsb|lm` (`refine3D_local refine=prob_neigh`, one state, euclid,
    trs > 0; nspace 1,000 unless given, nspace_sub = nspace). With it on, the master writes the B-spline references
    beside the polar model (`cart_refs_needed`), each part reads them into a calculator of its own and builds
    `build%pcont` on it and `pftc`; `strategy3D_srch` gives the committed pose (not the extras) the 5-variable solve of
    `pftc_shsrch_grad`'s view mode (`'lbfgsb'` or `'lm'`, tilts bounded by 2 sqrt(4 pi/nspace) rad, the joint mode's shift
    prior, commit only on improvement over its own B-spline re-score); the matcher applies `commit_view` after the
    search (R = euler2m(assigned)·R_refᵀ·T(w)·R_ref, `proj` the nearest grid direction), takes such a particle's sigma2
    from the final evaluation (`calc_sigma2_view`) and writes the BENCH counts (`match3D view solves`, `sections cut`,
    `cost evaluations`, per solve).
  - **Exit 1, the no-change case and `fills` pass.** `cont_view=no`, Release, `nthr=1`: the 63 identity outputs
    byte-identical to Phase 0's base run, every other file in the declared exclusions (12 NORM_SAME;
    `scratch/logs/cmp_p2_prob_neigh_final.txt`, 260.5 s; also on the first Phase 2 build, `cmp_p2_prob_neigh_1.txt`).
    `fills` (Release, 8 threads): 48 rows, no differing entry, the table identical to the base's but for wall time
    (`scratch/logs/fills_p2_release.tsv`, 307.8 s).
  - **Exit 2, the bank and nspace_sub.** With the view on the bank holds exactly nspace = 1,000 reprojections per half:
    `reprojection_model_{even,odd}.bin` = 16 + 1,000 × pftsz × nk × 8 bytes (66,816,016 in the `lbfgsb` smoke run, the
    base's per-section size at 20,000 directions times 1,000), `NSPACE PROJECTION DIRECTIONS | 1000` in `run.log`,
    `match3D nspace : 1000` in BENCH. An explicit different nspace_sub is refused by the commander
    (`scratch/logs/refusals.txt`, case 2) and by the parameters (`test_cont_view`).
  - **Exit 3, the recovery** (`test_view_recovery`, `simple_strategy3D_inplane_tester.f90`; the 4-blob phantom, lp 6 Å
    fixed before the first run, cell 500 of the 1,000-direction spiral, a noiseless KB projection 3° off it (tilt
    direction 60°, the first 30° step keeping the cell nearest), 17.3° in-plane, (1.4, −0.6) px; seeded at the cell, its
    nearest grid rotation, zero shift; driven through `refine_assignment_continuously`, `assign_ori` and `commit_view`):
    view error `lbfgsb` 0.0545° (Debug 0.0596°), `lm` 0.0596° (Debug 0.0596°) against the floor 0.25° (the plan's,
    guessed); with the view off it keeps 3.0000° (asserted within 0.1°); `proj` is the nearest grid direction to the
    committed pose. Sections per solve on the fixture: `lbfgsb` 65 (67 evaluations; Debug 39/40), `lm` 4 (4; Debug 5/5)
    (`scratch/runs/p2_final/unit_pftc_align2D3D_{release,debug}.log`).
  - **Exit 4, sigma2.** The residual of the final evaluation (`polarft_cont_calc%sigma_contrib`) equals an independent
    recomputation (B-spline section at the solution's tilts, shifted, rotated, minus the particle) within 5e-16 of a
    2.7e-10 maximum (1e-4 asserted), and is the cost's residual: off the solution on a grid rotation the loss equals
    nrots × Σ_k w_k residual_k over the particle's weighted power, 1.90147e-2 against 1.90148e-2 (1e-4 asserted).
  - **Exit 5, refused off-route** (`scratch/logs/refusals.txt`, `simple_exec` on the fixture project, each stopping before
    any worker): `refine3D_local refine=cont`, nspace_sub ≠ nspace, `cont_view=maybe`, `nstates=2`, `refine3D_auto`. The
    parameters' route checks in forked children (`test_cont_view`, `unit_core` 1138 assertions): a bad value,
    nspace_sub ≠ nspace, `refine=cont`, two states and `objfun=cc` stop; `lm` and `lbfgsb` on the route construct; another
    program ignores the key. Plain `refine3D refine=prob_neigh` accepts it (deviation, decision 1 below).
  - **Exit 6, sections per solve on the fixture:** Exit 3.
  - **Exit 7, smoke runs** (the identity command at `nthr=8`, nspace = nspace_sub = 1,000, `scratch/runs/smoke_*`,
    repeated on the final tree in `scratch/runs/smoke2_*` with the same counts and resolutions): all complete.

      | `cont_view` | wall s (run 1 / repeat) | FSC 0.143 / 0.5 (Å) | sections per solve, iterations 1-4 | evaluations per solve, iteration 1 |
      | --- | --- | --- | --- | --- |
      | `no` | 30.2 / 29.5 | 3.47 / 3.98 | - | - |
      | `lbfgsb` | 53.0 / 53.2 | 3.17 / 3.60 | 53.4, 55.6, 55.3, 54.4 | 56.4 |
      | `lm` | 35.1 / 34.8 | 3.17 / 3.60 | 15.3, 9.8, 10.0, 9.7 | 14.7 |

      Iteration 1 solves every particle (1,000), the others the 500 sampled. Load 2.3-9.2, my sequential chain and the OS
      daemons: the times are secondary. A smoke run, not a verdict; the view-off arm at the same nspace is C's shape.
  - **Gates.** Fast gate 15/15: Release 17.1 s, Debug 53.0 s (`scratch/logs/fastgate_p2_{release,debug}.log`). Suites,
    Release and Debug: `unit_pftc_align2D3D` 6471 (Phase 1: 6457; the recovery test 14), `unit_core` 1138,
    `unit_cart_align3D` 1178 (`scratch/runs/p2_final/`). T4 and its negative controls pass.
  - **Decisions under delegation.**
    1. Plain `refine3D` accepts the key on `refine=prob_neigh` with the same checks: it is the stage program
       `refine3D_local` runs in process and its workers run as, and only the router (`exec/simple_exec_refine3D.f90`,
       not in the table) can tell a user from the stage (precedent: `cont_prior`). Every other program ignores it.
    2. nspace defaults to 1,000 under the key (decision 8's "(1,000)").
    3. The commit is applied by the matcher after the search (`assign_ori` and `strategy3D_prob` are not in the table):
       `strategy3D_srch` keeps the solution, `commit_view` overwrites the Euler triplet and `proj`.
    4. The B-spline references of a view pass go into a calculator owned by the pass, not `build%cftc`, whose references
       make each batch prepare Cartesian particle slots (`build_batch_particles3D`, not in the table): the first smoke
       runs stopped on "fewer particle slots than the batch".
    5. Only a committed (improving) view solve takes the final evaluation's residual; the others keep the bank's.
       `gen_sigma_contrib` was split into a loader and one body so the bank path keeps its operations.
    6. A non-improving view solve stores the B-spline re-score as `corr` (the contract's B-spline-only comparison).
    7. `cont_view` with `trs=0` is refused (it would only force the grid), after the independent review of the diff.
    8. The recovery fixture's tilt direction is the first 30° step that keeps the drawn cell nearest (precondition).
  - **Finding (base, not changed).** `calc_sigma2_pftc` applies the matcher's stored-frame shift increment as the native
    shift of `gen_sigma_contrib`: with the view off on the recovery fixture, the residual so computed is 1.17e-2 against
    the joint cost's 6.26e-3, the native shift giving 6.26e-3 (`scratch/runs/p2_tests/rec5.log`, the diagnostic line).
    It inflates every polar refinement's sigma2 at a nonzero in-plane angle; view-committed particles take the cost's
    residual instead, a confound for the PorF A/C against H comparison. Reported for a ruling (report, section 7).
  - What it means, in short: `cont_view=lbfgsb|lm` refines the drawn pose's view on B-spline sections with the polar
    search on a 1,000-direction bank, recovers a 3° view error to 0.06° on a noiseless fixture, runs end to end with LM
    at about a quarter of L-BFGS-B's sections per solve, and leaves `cont_view=no` byte-identical; the PorF read is
    Hans's.
- 2026-11-19, session (FWM-S116118, on `f0211f929`): the run's `review/full.diff` applied to the checkout (one
  conflict, `simple_parameters_tester.f90`: both new tests kept). Hans's rulings on the report: the sigma2 frame
  is corrected; plain `refine3D` accepting the key stays; the deviations stand; the skill text of report section 6
  is applied; Hans runs Phase 3.
  - **Sigma2 frame fixed.** `calc_sigma2_pftc` rotates the stored-frame increment into the particle frame at the
    evaluated rotation before `gen_sigma_contrib`. It changes every polar refinement's sigma2 (2D and 3D) at a
    nonzero in-plane angle and shift; the no-change case of this plan no longer holds against the base for that
    reason alone. Tests: `sigma2_anchor` `test_shift_frame` (the residual is taken at the particle-frame shift;
    the fixture discriminates) and the recovery's view-off branch, whose diagnostic is now an assertion (bank
    residual at the native shift = joint cost within 1e-4; the unrotated increment differs).
- 2026-11-20, session (on `4ac019e30`; Mazharul Hossain's review, Hans's ruling: the complete fix). When a refined
  extra in the committed pose's basin beats it, `resolve_extra_poses_view` replaces the committed pose, but the
  dynamic centring (`record_shift_state`) still took the original candidate's shift increment and its cell's
  posterior width. A move of the replacing extra that reached the search bound then failed to widen the next
  search. The insertion table's extras now carry their candidate's posterior width (`prob_insert_pose%sdev`, from
  `prob_candidate%sdev` in `insert_table_add`), and the replacement hands its refined increment and that width to
  `record_shift_state`. Only `cont_view` with `topk > 1` and dynamic centring changes. Test:
  `test_neigh_pooled_extras` checks every extra's width against its candidate's.
