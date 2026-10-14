# The probabilistic search: paying less for the same table

Status: refactoring plan, 2026-10-13. No code yet. It follows the ribosome
tally of `doc/implementation_notes/planned/particle_caching_3D.md`, which
put the probabilistic pre-step at 46% of a `solve3D` chunk run, and the
code review of the same day (folded into section 2 below). The approach is
unchanged: the same exhaustive table, the same likelihood sampling, the same
assignment. What changes is where the shift ramp is applied, how the
sub-pixel shift is found, and how the table loop walks memory.

Rulings sought (section 9): none given yet. Phase 0 (section 10) runs
without any.

## 1. Aim

- Cut the time of the `prob_tab` parts without changing which references,
  rotations and shifts the table records, except where section 3.2 says so
  and the validation of section 6 accepts it.
- Give the pre-step benchmark files of its own, so its time is split into
  named parts instead of being one number in the refine3D strategy file.
- Decide the table-loop change by a measurement, not by an argument.

## 2. What the pre-step costs today

Measured (the ribosome chunk, consolidated layout, 29,101 particles, three
states, `nspace` 1,000, 5 parts of 16 threads, box 336): the pre-step is
3,033 of 6,572 s over 78 iterations, 37 s per iteration at `box_crop` 132,
50 s at 160 and 138 s of a 190 s iteration at 182. Inside it only the whole
is timed (`rt_prob` in `simple_refine3D_strategy.f90`).

Counted from the code, per sampled particle, with more than one state, at
`box_crop` 182 (`nk` 60 shells, `pftsz` 254, `nrots` 508):

- **The table** (`fill_tab_range`, `simple_eul_prob_tab.f90`;
  `process_particle_with_shift`). One thread per particle. The shift seed is
  single-state only (`l_sh_first`), so every one of the 3,000 references is
  scored at zero shift by `gen_prob_likelihood_euclid_val`: the memoised
  contraction `ft_ctf2 · ft_ref2 − 2 ft_ptcl_ctf · conj(ft_ref)` over
  `nk × (pftsz+1)`, one real inverse FFT of length `nrots`, and a bounded
  heap for the best `nsample` rotations (`sample_likelihood_dist`, already
  O(n log nsample)). Call the contraction plus inverse FFT the base kernel:
  about 0.32 Mflop, reading 248 KB of reference and 248 KB of particle
  memos. 3,000 base kernels per particle.
- **The refinement** (`refine_best_refs`). The `projs_ns` best references,
  `projs_ns = athres/180 × nspace` with `athres` the smaller of
  `prob_athres` (10°) and the mean angular move of the sampled particles:
  17 to 55 references. Each gets `grad_shsrch_minimize`
  (`simple_pftc_shsrch_grad.f90`) in alternating mode with `coarse_init`
  on, which the `xy_in` path of the single-state case turns off. Per
  refined reference:
  - one base kernel at zero shift;
  - the coarse grid: `coarse_num_steps`² = 25 calls of `gen_objfun_vals`,
    each shifting the *reference* (`gen_shmat4aln`, then `nk` complex FFTs
    of length `nrots` on `plan_fwd1_many`) before the contraction: about
    1.7 Mflop, five base kernels, the FFTs four fifths of it;
  - up to `max_evals` = 5 alternations of L-BFGS-B (`maxits_sh` = 60) and a
    shifted `gen_objfun_vals` each; every L-BFGS-B evaluation
    (`gen_corr_grad_for_rot_8` → `gen_euclid_grad_for_rot_8`) copies the
    reference to double precision (244 KB) and sums residual and gradient
    directly over `nk × pftsz` with the ramp regenerated: about 0.6 Mflop in
    double. The evaluations per run, E, are counted by the search object
    (`profile_objective_evals`, `profile_gradient_evals`) and printed
    nowhere.

  In base kernels per refined reference: 1 + 25 × 5.3 + alternations ×
  (2E + 5.3), about 250 to 450 for two to three alternations and E between
  20 and 50. Per particle: 4,000 to 25,000, against 3,000 for the table.

The consistency check: 5,000 particles × about 13,000 base kernels × 0.32
Mflop is 2 × 10¹³ flop per iteration, 1.9 Gflop/s per thread over 138 s on
80 threads, a plausible rate for these loops. The counts and the measured
time agree, so the pre-step is these two loops plus, of the order of 10 to
15 s per iteration, the second process launch, the second reprojection-model
read and memoisation per part, the table write and read, and `ref_assign`.

The same grid runs once per particle in the matcher
(`grad_shsrch_first_obj` in `simple_strategy3D_srch.f90`, then the peaks
with `xy_in`), in the 2D search (`simple_strategy2D_srch.f90`), in the 2D
table (`simple_eul_prob_tab2D.f90`) and in the neighbourhood table
(`simple_eul_prob_tab_neigh.f90`). The multi-state 3D table is the one place
that runs it per refined reference.

## 3. Target design

### 3.1 The coarse grid on the particle side (exact)

For a unit-modulus shift ramp S and a real CTF,
`|X·S(s) − CTF·R|² = |X − CTF·R·S(−s)|²` at every polar sample, so
shifting the particle by s gives the same Euclidean distances, at every
rotation, as shifting the reference by −s; `ft_ctf2` and the particle's
weighted power are unchanged (|S| = 1). The 25 grid shifts are fixed by
`trs` and the same for every reference of a particle.

- The search object pre-shifts the particle at the 25 grid positions once
  per particle: per position a ramp and the `nk` FFTs over rotation of
  `X·S·CTF`, kept as 25 memoised `ft_ptcl_ctf` sets in a per-thread
  workspace of the polar calculator (25 × `(pftsz+1) × nk` complex, 3 MB
  per thread at `box_crop` 182). They are built when the object is set to a
  new particle (`set_indices` with a changed particle) and reused for every
  reference refined for that particle.
- `coarse_search_opt_angle` and `coarse_search` evaluate each grid point
  with a variant of the unshifted branch of `gen_raw_euclid_vals_impl` (and
  of the CC path in `gen_corrs`) that takes the pre-shifted particle memo
  in place of `ft_ptcl_ctf(:,:,i)`: a base kernel per grid point.
- The grid's shift is reported as today, as a reference shift: the particle
  memo at grid position s stands for the reference shift −s. The sign lives
  in one place, the workspace builder.
- Per refined reference the grid costs 25 base kernels instead of 25 × 5.3;
  the 25 ramps and FFT sets are paid once per particle. For `projs_ns` = 30
  that removes about 1,000 of the roughly 3,400 Mflop per particle of
  section 2.
- Every consumer of `coarse_init` gains; the per-particle case (matcher,
  2D) gains less, since there the grid was already once per particle.

### 3.2 The sub-pixel step from the grid (approximate, validated)

The grid spacing is `2·trs/5`. Around the grid minimum the cost is smooth in
the shift, and the 3 × 3 neighbouring grid values are already computed: a
quadratic fit through them gives the minimum to a fraction of the spacing
in one step. The values used are those at the rotation that is best at the
centre point (the kernel returns all rotations for every grid point, so
this costs nothing). One exact kernel at the fitted shift then gives the
distance that the table records; the fit is never trusted for the value.

What follows the fit is a ruling (section 9): either L-BFGS-B starts from
the fitted shift and converges in a few evaluations, or the table fill
records the fitted shift and leaves the continuous refinement to the
matcher, which refines the assigned candidate anyway. Either way this
replaces alternations × 2E ≈ 100 to 300 base kernels per refined reference
with about one, or with a few.

This changes the recorded shifts by less than the grid spacing, so it is
validated as section 6 says, not assumed.

### 3.3 The table loop and memory (decided by Phase 0)

The table loop streams, per particle and thread, the whole memoised
reference bank of the particle's parity (3,000 × 248 KB = 744 MB) at under
one flop per byte; the particle's own memos stay in cache. Two forms:

- **Blocking.** A thread takes a block of B particles and, for each
  reference, evaluates the B pairs before moving to the next reference, so a
  reference's memos are loaded once per block. The per-pair work (contraction,
  inverse FFT, heap, candidate record) is unchanged; only the loop order in
  `process_particle_with_shift`/`process_particle_without_shift` changes.
- **The product form.** For each of the `pftsz+1` rotation frequencies the
  contraction is a complex matrix product over the `nk` shells (particles ×
  shells times shells × references); a BLAS routine runs it at a multiple of
  the hand loop's rate. More invasive: a new kernel in the polar calculator
  with a block of particles as input.

Phase 0's thread-scaling run decides: if one `prob_tab` part scales linearly
from 1 to 16 threads the loop is compute-bound and only the product form
helps; if the scaling flattens it is bandwidth-bound and blocking alone
helps. Blocking first in either case; the product form only if blocking
leaves the table the larger part.

### 3.4 Benchmark files

- `PROB_TAB_BENCH_ITER###_PART##.txt`, written by `exec_prob_tab` in the
  matcher's format (`CONTEXT`, `TIMINGS`, `% accounted for`): project read,
  reprojection-model read, reference memoisation, particle preparation, the
  score loop, the refinement, the table write; and in the context the
  counts: references, sampled particles, `projs_ns`, the search object's
  objective and gradient evaluations summed over threads, the number of
  alternations, the grid evaluations.
- `PROB_ALIGN_BENCH_ITER###.txt`, written by `exec_prob_align`: sampling and
  project write, scheduling (launch to last part finished), table reads,
  `ref_assign`, assignment write.
- The score loop and the refinement are timed inside the parallel loop with
  per-thread accumulators summed after it; the counters come from the
  per-thread search objects after the loop. Names in
  `simple_refine3D_fnames.f90` next to `refine3D_bench_fname`.

## 4. What does not change

- The table's contents: every reference scored, the same rotations kept by
  the same heap, the same likelihood sampling, the same assignment.
- The single-state seed path and the matcher's seed-then-refine pattern.
- The reference-shift convention in every recorded shift.
- The 2D table and the neighbourhood table, beyond gaining 3.1 through the
  search object.

## 5. Order of work

1. Benchmark files (3.4), the thread-scaling run, the baseline: Phase 0.
2. The particle-side grid (3.1), exact; tested for equality against the
   reference-side grid: Phase 1.
3. The sub-pixel step (3.2), after its ruling, validated against ground
   truth: Phase 2.
4. The table loop (3.3), if Phase 0 says it is bandwidth-bound or if the
   table is the larger part after Phase 2: Phase 3.
5. Documents and close: Phase 4.

Each step is measured on the same case before and after; a step is kept if
it is faster beyond the run-to-run spread and the acceptance holds.

## 6. Tests

- **The particle-side grid.** In `simple_pftc_inplane_tester.f90`: for
  random particle and reference polar transforms with CTF, the distances
  over all rotations from the reference-side evaluation at shift s and from
  the particle-side evaluation at −s agree to single-precision rounding
  (relative 1e-5), under both objectives; and the grid minimum, its rotation
  and its shift from `coarse_search_opt_angle` are identical before and
  after the change on the same inputs (same evaluation order, same
  arithmetic up to the rounding above; the test asserts the same grid point
  and rotation).
- **The sub-pixel step.** A tester case on a synthetic particle with a known
  shift: the fitted shift lies within 0.2 pixel of the L-BFGS-B result from
  the same grid point, over a range of true shifts and noise levels; and
  the workflow validation of the test protocol: on a small noisy simulated
  set (2,000 particles or fewer) the orientation and shift errors against
  ground truth before and after lie within the baseline spread; on
  beta-galactosidase the final FSC 0.143 resolution of every state lies
  within the larger of the baseline spread and 3% of the baseline.
- **The table loop.** The per-pair minimum distance and its rotation are
  identical before and after blocking (the deterministic part of
  `sample_likelihood_dist`); the sampled rotation draws `ran3` in a
  different order per thread, so the sampled picks are compared as
  distributions (fraction of particles whose assignment changes, against the
  run-to-run fraction).
- **Benchmark files.** Their `% accounted for` is above 90% on the timed
  cases.
- `python3 scripts/check_test_registry.py .` whenever tests or their
  registration change.

## 7. Documents

- `doc/policies/3D/refine3D_policy.md`: the pre-step's benchmark files; the
  particle-side grid stated as an implementation of the same objective; the
  sub-pixel step if ruled.
- `.github/skills/simple-main-pftc/SKILL.md`: the pre-shifted particle
  workspace as part of the memoisation design (text proposed in the report,
  applied on approval).
- `doc/implementation_notes/planned/particle_caching_3D.md`: a pointer to
  this plan replaces its "what the measurement puts first" paragraph on the
  pre-step.
- On completion this plan moves to `doc/refactoring_notes/completed/` with
  its report.

## 8. Consequences

- Memory: 3 MB per thread of pre-shifted particle memos; nothing else.
- The prob-stage iteration on the ribosome chunk: by the counts of section
  2, the grid change alone is worth a quarter to a third of the pre-step and
  the sub-pixel step another quarter to a half; the table-loop change is
  worth something only if Phase 0 shows the loop bandwidth-bound. These are
  estimates; Phase 0 replaces them with numbers before anything is built.
- The matcher's discrete search and the 2D paths gain the grid change for
  free, per particle rather than per reference, so less.

## 9. Open rulings

1. After the quadratic fit: L-BFGS-B from the fitted shift, or no
   continuous refinement in the table fill at all (the matcher refines the
   assigned candidate).
2. Per-state seeds for more than one state: the refine3D policy forbids a
   seed from the particle's previous *state*; a seed estimated per state at
   the particle's previous projection direction would remove the grid from
   the multi-state table as the single-state path already does, at one grid
   per state per particle. This changes the search and is not in the plan
   until ruled.
3. The search object's evaluation counters: printed in the benchmark files
   (this plan) rather than deleted as unused.
4. Deferred until the global sigma2 refactoring has landed: with one shell
   weight for all particles, the `ft_ctf2 · ft_ref2` half of the contraction
   depends on the particle only through its CTF, so for the particles of one
   micrograph it is one computation per reference; worth about a third of
   the table's contraction on the ribosome chunk (about four sampled
   particles per micrograph). Not in this plan.

## 10. Execution

Carried out as an unattended run on the Oracle Linux test machine (arun
run `prob_search_speedup`) for the code, the tests and the relative
timings; the ribosome chunk timing on the production server is the
maintainer's, before Phase 0 closes and after Phase 2. Nothing is committed;
the maintainer reviews the complete diff. Every phase records itself in
section 11.

Common to every phase:

- **Data.** Beta-galactosidase (5,513 particles, box 256 at 1.275 Å per
  pixel, `mskdiam=180 pgrp=d2`, 50 selected class averages), copied fresh
  for every run; the originals are never written.
- **Execution settings.** `nparts=4 nthr=6` distributed, identical in every
  phase; timed runs on an out-of-tree Release build with tests built,
  nothing else heavy running.
- **Cases.**
  - **S.** `solve3D_cavgs`, single state; kept as the common start.
  - **P1.** `solve3D nsample=2500 cavg_ini_ext=yes` from S, single state:
    the seed path.
  - **P3.** The same with `nstates=3`: the multi-state table, the case of
    the ribosome run.
  - **T.** One `prob_tab` part of a P3 prob-stage iteration, rerun alone on
    the same particles at `nthr` 1, 2, 4, 8 and 16 (the saved project and
    reprojection model of that iteration; `part=1`).
  - **W.** The simulated workflow tests with the Euclidean objective, for
    the ground-truth errors.
- **Acceptance of a repeated case.** As the test protocol: P1 and P3 by the
  final FSC 0.143 resolution of every state within the larger of the Phase 0
  spread and 3%; W by orientation and shift errors within the Phase 0
  values plus the same margin.
- **Gates.** The fast gate (`./compile_debug.sh`), the CTest unit suites of
  the touched areas, and the registry check whenever tests change.

### 10.1 Phase 0: benchmark files, baseline, thread scaling

1. Build the tree as it is (debug with the fast gate, Release for timing).
2. The benchmark files of 3.4, the only source change of this phase, in
   SIMPLE's benchmark conventions (`L_BENCH_GLOB`, `tic`/`toc`, no new test
   arrays).
3. P1 and P3 twice each (the spread), W once; record per case the command
   line, the wall time, the final resolutions, and from the new files the
   score-loop and refinement times, E, the alternations and the grid
   evaluations per prob-stage iteration.
4. T at the five thread counts; record the score-loop and refinement time
   at each.

Exit: the files appear with `% accounted for` above 90%; the split and the
scaling are in section 11 with the decision for Phase 3 (compute- or
bandwidth-bound) written down.

### 10.2 Phase 1: the particle-side grid

3.1 in the search object and the polar calculator; the tester case of
section 6. P3 repeated: the acceptance holds and the refinement time is
compared with Phase 0. P1 and the matcher's cases are not expected to
change beyond the spread; P1 repeated once to confirm.

Exit: the equality test passes; P3's before-and-after is in section 11.

### 10.3 Phase 2: the sub-pixel step

After ruling 1. 3.2 in the search object; the tester case and the
validation of section 6 (W, then P3 and P1). Kept only if the acceptance
holds.

Exit: W within acceptance; P3 and P1 within acceptance; the
before-and-after in section 11.

### 10.4 Phase 3: the table loop

If Phase 0 showed the loop bandwidth-bound, or if after Phase 2 the score
loop is the larger part: blocking in `fill_tab_range`, B chosen from T at
8, 16 and 32; the test of section 6. The product form only if blocking
leaves the table the larger part, as its own entry.

Exit: the test passes; T and P3 before-and-after in section 11.

### 10.5 Phase 4: documents and close

The documents of section 7 except the skill file and the generated indexes;
P1 and P3 once more on the final build; the plan moved to
`doc/refactoring_notes/completed/` with the report beside it, which reads on
its own: what changed and why, the evidence against Phase 0 in tables, every
decision taken and where it is recorded, the proposed skill text, what is
left open.

### File table

| File | Change | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_prob.f90` | Benchmark files of `exec_prob_tab` and `exec_prob_align`. | 0 |
| `src/defs/simple_refine3D_fnames.f90` | `prob_tab_bench_fname`, `prob_align_bench_fname`. | 0 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | Per-thread timers and counters in `fill_tab_range` with getters; blocked loop order. | 0 3 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | Counter getters; pre-shifted particle workspace on `set_indices`; grid on the particle side; the quadratic step. | 0 1 2 |
| `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_memo.f90`, `simple_polarft_corr.f90` | Per-thread pre-shifted particle memos; the evaluation variants taking a particle memo (Euclid and CC); a blocked table kernel if the product form is built. | 1 3 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | The equality and sub-pixel cases. | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | Nothing unless the blocked loop needs a block-wise candidate record. | 3 |
| `doc/policies/3D/refine3D_policy.md`, `doc/implementation_notes/planned/particle_caching_3D.md` | Section 7. | 4 |
| `doc/refactoring_notes/completed/prob_search_speedup.md`, `prob_search_speedup_report.md` | The plan moved on completion, and the report. | 4 |

## 11. Progress

No phase has run yet.

## 12. Code map

| Area | Files |
| --- | --- |
| Pre-step orchestration | `src/main/strategies/parallelization/simple_refine3D_strategy.f90` (`distr_execute_iteration`, `rt_prob`), `src/main/commanders/simple/simple_commanders_prob.f90` (`exec_prob_align`, `exec_prob_tab`) |
| The table | `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` (`fill_tab_range`, `refine_best_refs`, `ref_assign`), `simple_eul_prob_tab_utils.f90` (`sample_likelihood_dist`, `calc_num2sample`) |
| The shift search | `src/main/pftc/simple_pftc_shsrch_grad.f90` (`grad_shsrch_minimize`, `coarse_search_opt_angle`, `coarse_search`, the counters) |
| The kernels | `src/main/pftc/simple_polarft_corr.f90` (`gen_raw_euclid_vals_impl`, `gen_corrs`, `gen_prob_likelihood_euclid_val`, `gen_euclid_grad_for_rot_8`), `simple_polarft_memo.f90` (`memoize_refs`, `memoize_ptcls`) |
| Other grid users | `src/main/strategies/search/simple_strategy3D_srch.f90`, `simple_strategy2D_srch.f90`, `probabilistic/simple_eul_prob_tab2D.f90`, `simple_eul_prob_tab_neigh.f90` |
