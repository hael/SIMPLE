# In-plane quadrature insertion in 2D (`topk_inpl`)

Status, 2026-11-13: implemented and gated, uncommitted, awaiting Hans's validation of the effect; the report is
`inpl_quadrature_insertion_2D_report.md` beside this note. Revised the same day after review R1 (section 10).
Base: `db20cde75`.
Supersedes: none. Depends on: the 2D batch pipeline (`eul_prob_tab2D`, validated: yes, in production) and
the pattern of the 3D weighted insertion's side channel (`prob_insert_poses`, validated: yes, production
since 2026-10-30; its selector and cone are not reused). Rulings 1 to 4 and 14 are Hans's (2026-11-13);
the rest are proposed. The plan's gates are the byte identity of `topk_inpl=1` with the current code and
the unit assertions of the enabled path; the larger-scale testing of the effect is Hans's, after delivery
(ruling 4).

## 1. Why

Top-view class averages of small membrane proteins (SLC transporters, hemoglobin-sized particles) from
`solve2D` show spurious radial artifacts ("spokes"). They persist with the gold standard from stage 2, in
runs with a limited `lpstop`, and independently of the annealing tail (Hans, 2026-11-13); Joe Caesar has
seen them appear as early as stage 2 and never leave. The sigma annealing proposal of 2026-11-10 (not in
the repository) does not address them directly: every noise model in it is a function of radius alone,
hence circularly symmetric in Fourier space, and cannot produce angular structure.

The mechanism this plan acts on: a top view in its micelle is close to rotationally symmetric, so the
in-plane angle of a particle carries almost no signal and is decided by the particle's noise against the
reference's asymmetric residue. The particle's noise then enters the half's sum at that angle and the residue
grows; a loop whose gain exceeds one keeps what it has made, which is Joe's observation. The gold standard
keeps the half-FRC honest but does nothing within a half. Per half and per shell the restored averages do
not depend on sigma2 given the poses (section 2, with its assumptions), so the noise model enters this loop
through the search: which rings vote on the angle, and how sharply.

**What the insertion is hypothesised to do.** The particle enters the sums at a quadrature of its in-plane
posterior, conditional on its committed class and its polished shift, instead of at one drawn angle. The
hypothesis is that this lowers the fresh asymmetric noise the loop is seeded with: a `K`-point quadrature
of a flat row, at equally spaced angles, annihilates every angular harmonic of the particle's noise except the
multiples of `K` (review R1, point 1); with `K` above the band's angular Nyquist at the rim (ruling 15) no
resolved harmonic is one of those multiples. The points are snapped to the rotation grid, which leaves a
residual unless `K` divides the number of rotations (review R2, point 2). It is not a rotational average over the full circle, and it does not integrate class or
shift uncertainty. Whether the enabled update amplifies or attenuates a small asymmetric reference
perturbation is a separate question: the quadrature is taken at the particle's own temperature
(ruling 14) while the search draws at `prob_beta` times it, and the response of a posterior expectation to
a small perturbation scales with the temperature (R1, point 4). The plan does not claim the loop's gain is
unchanged; the perturbation measurement of section 7 is the way to settle it, if Hans wants it run. The
intended benefit must survive the terminal hard pass through better references and assignments: smoother
stage averages alone prove nothing (R1, closing).

## 2. What the code check changed

Checked at `db20cde75` with `rg`; these premises are fixed by the code.

- **The in-plane draw is already hot.** The fill draws each class's in-plane index at
  `ptcl_posterior_temp`, which is `prob_beta` times the particle's temperature
  (`simple_polarft_corr.f90`, `ptcl_posterior_temp`), and the class draw is at the same factor
  (`assign_rows_dense`, `assign_rows_sparse`). The 3D extras are weighed at the particle's own temperature,
  untempered (`select_extras`, `temperature_or_argmin(self%temps(i_loc))`); the 2D quadrature does the
  same (ruling 14).
- **The 3D extras average over direction within a cone**, each entry with its own in-plane index. A 2D
  mirror over the row (classes) would be soft class membership and would not touch the angle; a cone in
  angle (rim displacement below the resolution) would exclude the whole undecided circle of a top view.
  Hence ruling 5: the angle, the whole row, no cone.
- **The committed angle is continuous.** After the joint polish `store_continuous_e3` keeps a fractional
  angle and `assign_ori` writes it when the class and grid index are the committed ones
  (`strategy2D_srch.f90`, `has_continuous_e3`); otherwise e3 is `360 - get_rot(irot)`. A quadrature on
  the grid alone would replace a polished angle by its bin centre (R1, point 2): ruling 16.
- **The row at the committed pose has a public evaluator.** `gen_euclid_dists_bank(iref, iptcl, bank,
  dists)` returns the `nrots` losses of the particle shifted once in its own frame
  (`memoize_shifted_ptcl`) against the reference; the crvec helpers it wraps are private to the calculator
  (R1, implementation note). Ruling 7 uses it.
- **The 2D row holds one in-plane index per class** (`prob_dense_entry%inpl`); the full rotation row is
  computed in the fill (`gen_posterior_val`, `draw_posterior_crvec`) and not kept, so it is evaluated again
  at the committed pose (ruling 7).
- **The sums rotate at insertion with an optional weight.** `gen_fplane4rec` makes one unrotated plane per
  particle with shift and CTF applied; `stack_accumulate_fplane(e3, fpl, icls, weight)` rotates and grids
  it, the data term and the CTF-squared term alike. An extra angle is one more gridding of the same plane:
  no second read, FFT or phase ramp.
- **Sigma2 and the restored averages, with the assumptions.** The plane carries the data term and the
  CTF-squared term whitened by the same curve (`simple_image_ctf.f90`, `gen_fplane4rec`: "CTF*y/sigma2 and
  CTF^2/sigma2"). Within one half, at one shell, with sums of one iteration, the curve cancels in the
  ratio, and the Wiener step (`add_invnoisepower2rho`) forms its terms from the CTF sums and the FRC. The
  cancellation is not exact beyond that (R1, point 5): the KB window mixes neighbouring radii, the merged
  average adds the even and odd sums before normalisation (`cavger_restore_cavgs`, "even + odd"), so it
  weights the halves by their curves, and the carry-over blend (`cavg_sums%blend`) adds sums whitened under
  different generations. The algorithm doc states the cancellation under these assumptions and nothing more.
- **The row's shift frame.** `gen_euclid_crvec(iref, iptcl, shift)` applies the shift to the reference
  before the rotation, so across rotations it describes a particle shift that rotates with the reference.
  The bank form shifts the particle once in its own frame and gives a row consistent across the circle; the
  sums apply the stored shift to the particle before rotation (`cavger_update_sums`, `precs(i)%shift`), the
  same frame.
- **Where the search commits.** `strategy2D_prob%srch` takes the table's class, in-plane index and shift,
  runs `refine_selected_continuously`, stores and assigns; nothing after it touches the pose. The quadrature
  goes there, after the polish (ruling 7). `strategy3D_prob%resolve_extra_poses` is the pattern.
- **Which stages have a table.** Stage 1 runs `snhc_smpl` (classic, no table); stage 2 `prob`; stages 3
  to 5 `prob_neigh`; stage 6 `prob` (`set_refine2D_stage_search_policy`); the terminal pass is
  `refine=greedy` (`execute_terminal_pass`). The insertion exists in stages 2 to 6 only.
- **Keyed draws.** Every draw of the probabilistic routes is keyed on the run's scope and the particle
  (`rnd_keyed(rnd_draw_key(scope, iptcl, family))`, families `RND_FAM_*` in `simple_rnd.f90`); a new draw
  takes a new family, never `ran3`. Ruling 15's offset is such a draw.
- **No polar class accumulator exists**; the sums are Cartesian.
- **The readout's `n_eff` is over classes** (`neff_readout_offer` on the class distances), not angles.
- **`topk` is 3D only** in the interface (`simple_ui_solve3D.f90`, `simple_ui_refine3D.f90`); the 2D
  programs are all in `simple_ui_refine2D.f90` (refine2D, solve2D, solve2D_chunks). A chunk's solve2D command line
  is built from scratch with an explicit key list (`init_chunk_clustering`), so a solve2D key reaches the chunks
  only where that list forwards it (review R2, point 1; the first version of this note assumed otherwise).
- **The 2D bench** times the sums separately (`rt_cavg_interp_splat` in `restore_class_averages_for_batch`),
  so the cost of the extras is visible where it lands.

## 3. Contract

- **Off state.** `topk_inpl=1` (given explicitly; the default is 64, ruling 17) is today's run bit for bit at one thread: no row
  is evaluated, no side channel allocated, no readout field written other than zeros. This is the plan's
  gate and no ruling may waive it; the assertions of Phase 1b cover the enabled path.
- **The committed pose is the hypothesis.** For identical incoming references the search, its draws and
  keyed random streams, the polish, the project's class, e3 and shift, the sigma2 residual at the committed
  pose, the convergence statistics and every existing readout field are unchanged at any `topk_inpl`. The
  references change under `topk_inpl > 1`, so the run's trajectory changes from the second iteration on;
  "unchanged" is a statement about one pass, not about a run (R1, closing).
- **The concentrated-posterior limit reproduces today's insertion** (ruling 16): a particle whose
  quadrature collapses to one entry enters at its polished pose at weight 1, exactly as at `topk_inpl=1`.
  Enabling the feature never inserts a well-determined particle less accurately than today.
- **The sums stay thread-count independent.** A particle's entries are inserted in its slot of the bucket,
  in batch order (`cavger_update_sums`); weights of a particle sum to 1 and each weight applies to the data
  term and the CTF-squared term alike, so populations, the particle's mass and the carry-over's accounting
  are unchanged. The worker's contribution file format is unchanged.
- **Hard assignment at the end.** The terminal greedy pass stays (ruling 1): the final class averages are
  built at the committed poses and `make_cavgs` from the project reproduces them. The averages of the
  iterations of stages 2 to 6 under `topk_inpl > 1` are not reproducible from the project; the policy says
  so (ruling 1, Phase 1).
- **Scope of the device.** The weighting is internal (ruling 3): nothing downstream reads it, no file
  carries it, and 3D is untouched.

## 4. Decisions and rulings

1. (2026-11-13, Hans) The terminal greedy pass stays as it is. It is where the solution becomes the hard
   assignment. The solve2D policy states the difference in how the averages are generated: the stage
   iterations' averages take each particle at the quadrature of its in-plane posterior; the terminal pass's
   and `make_cavgs`'s take it at its committed pose.
2. (2026-11-13, Hans) The key is `topk_inpl`. `topk` keeps its 3D meaning and is not read in 2D.
3. (2026-11-13, Hans) The weighting is an internal trick against overfitting, not a model of the output.
4. (2026-11-13, Hans; restated the same day after R1) The plan carries no larger-scale testing: no fixture
   campaign, no real-data arms, no floors, no real-data read. The effect is tested by Hans after delivery.
   Unit assertions of the enabled path's numerical correctness (R1, point 3) are part of the plan
   (Phase 1b, section 7).
5. (2026-11-13, proposed; briefly replaced by ruling 11 the same day, restored with ruling 14) The extras
   are in-plane rotations of the committed class at the committed particle-frame shift, over the whole
   row. No cone in angle, no shift extras, no class extras.
6. (2026-11-13, proposed; its phase replaced by ruling 15, its temperature by 14, its collapsed case by
   16) Selection: a `K = topk_inpl` point quadrature of the row's posterior: weights
   `w_r = exp(-T (d_r - d_min))` over the `nrots` rotations, the cumulative sum `C`, the rotations
   containing `u_j = (j - 1 + v) / K`, `j = 1..K`, each at weight `1/K`, coincident rotations merged by
   summing. For `K > 1` the committed draw is not inserted separately (it is one sample of the same
   posterior; inserting it beside the quadrature would count it twice). Top-k is not used: on a near-flat
   row it returns the noise-best angles, the ones that make the spokes.
7. (2026-11-13, proposed; the evaluator corrected after R1) The row is evaluated in
   `strategy2D_prob%srch` after `refine_selected_continuously`, on the particle's thread:
   `memoize_shifted_ptcl(iptcl, best_shvec, bank)` then the public
   `gen_euclid_dists_bank(best_class, iptcl, bank, dists)`. An assertion of Phase 1b checks that the
   row's value at the committed discrete rotation agrees with `gen_corr_for_rot_8_shifted` at the same
   pose, with a non-zero particle-frame shift and an astigmatic CTF. The per-thread bank is allocated once
   per pass in `prep_strategy2D` (sized by the team), not per call.
8. (2026-11-13, proposed) The side channel is a new type `prob_insert_inpl` in `simple_eul_prob_tab_utils`
   (private components, `new`/`kill`, pattern of `prob_insert_poses`): per batch particle the merged
   angles as e3 (continuous, ruling 16) and their weights, laid out per batch, filled by the strategy, read
   by `cavger_update_sums` through a new optional argument. When a particle has entries the sum inserts
   them instead of the committed pose, each weight on both accumulators; without entries (stage 1, the
   terminal pass, `topk_inpl=1`, a particle without a valid row) it inserts the committed pose at weight 1
   as today. Readout: four fields, `inpl_nptcls` (particles with a quadrature), `inpl_npose` (entries after
   merging), `inpl_nsingle` (particles whose quadrature collapsed to the polished pose), `inpl_neff_sum`
   (sum over the particles with a quadrature of `exp(entropy)` of the row at the quadrature's
   temperature), written, read and added together (`write_prob_readout`, `read_prob_readout`,
   `add_prob_readout`), reported by the master as one `PROB_TAB_INPL` line: particles, entries per particle,
   the collapsed fraction, mean `n_eff`, and the stage's angular Nyquist `m_max` of ruling 15 beside `K`.
9. Replaced by ruling 14 (the temperature is the particle's own, not the draw's).
10. (2026-11-13, proposed) Default `topk_inpl=1` until Hans's testing decides; no upper bound, as `topk` has none
    (removed after the first delivery, which had a cap of 64 that nothing justified; Hans, 2026-11-13). The solve2D
    controller passes the user's value to stages 2 to 6 and
    deletes it on stage 1 and on the terminal pass; a direct `refine2D` caller gets it as given; the
    stream's pool runs the classic routes and never sees it.
11. Withdrawn the same day (a cone in angle, written on a misreading of "mirror 3D"; Hans: no cone).
12. Withdrawn with 11 (`top_posterior_extras` within the cone).
13. Withdrawn with 11 (the committed angle inserted beside the extras).
14. (2026-11-13, Hans) **The temperature.** The quadrature of ruling 6 is taken at the particle's own
    temperature, untempered (`temperature_or_argmin(posterior_temperature(wsqsum_ptcl_now, pftsz))`, not
    `ptcl_posterior_temp`), as the 3D extras are. This is the one element mirrored from 3D; the row is the
    whole circle (ruling 5) and the selection the quadrature (ruling 6). The search draws at `prob_beta`
    times that temperature; the inserted image is then the physical posterior's expectation, conditional on
    the class and the polished shift, approximated by `K` points, and the committed draw the hot hypothesis.
    The annealing tail leaves the quadrature's temperature where it is and cools the draws towards it.
15. (2026-11-13, proposed after R1 point 1) **The offset and the harmonic criterion.** The quadrature's
    phase `v` of ruling 6 is one keyed uniform variate per particle (`rnd_keyed(rnd_draw_key(scope, iptcl,
    RND_FAM_INPLQ))`, a new family), so a flat row gives `K` equally spaced angles at a phase that differs
    from particle to particle and the surviving harmonics (the multiples of `K`) carry random phases across
    the sum, instead of the same `K` grid angles for every particle. The angular-error criterion of the
    compressed insertion against the full posterior sum: for a flat row and equally spaced points the `K`-point
    rule equals the full sum on every harmonic that is not a multiple of `K`; the band's highest angular
    harmonic at the rim is `m_max = ceiling(pi * mskdiam / lp)` (`lp` the stage's matching limit in
    Angstrom, `mskdiam` the stage's mask diameter), so `K > m_max` keeps every resolved harmonic off the
    multiples of `K`. The points are snapped to the rotation grid (review R2, point 2), so the two agree
    exactly only where `K` divides the number of rotations (for a phase that is not a tie); otherwise harmonic
    `m` retains up to about `m 2 pi / nrots` (on 360 rotations `K = 64` retains 12.5 % of harmonic 40, `K = 32`
    13 % and `K = 16` 17 % of those up to 50). The readout prints `m_max` beside `K` per stage; Hans chooses
    `K` against it, knowing that `K > m_max` is necessary and not sufficient. For a peaked row the error of
    the midpoint rule on harmonic `m` is bounded by `m` times half the angular width of the widest stratum,
    which the stage's smallest `n_eff` makes small; this bound is stated, not gated.
16. (2026-11-13, proposed after R1 point 2) **The continuous-angle limit.** The quadrature's points are
    found on the grid; the point or points that fall in the committed bin are placed at the committed angle
    as the search commits it (`continuous_e3` when `has_continuous_e3`, else the bin's angle) and merged
    into one entry. A quadrature that collapses to one entry (every `u_j` in one bin, the concentrated
    posterior) inserts the polished pose at weight 1, bit for bit today's insertion for that particle; the
    readout counts these (`inpl_nsingle`). No continuous interpolation within other bins: a bin-uniform
    model of the posterior within a bin is wrong exactly where it would matter, at a peak.

Defaults taken without a question: the names `topk_inpl`, `prob_insert_inpl`, `PROB_TAB_INPL`,
`RND_FAM_INPLQ`.

17. (2026-11-13, Hans, after delivery and a run on beta-galactosidase with `topk_inpl=64`, class averages "almost
    perfect") **The default.** The parameter's default is 64 (Hans, the same day: one default every 2D method
    inherits, rather than commander defaults and callers pinned at 1), so refine2D, solve2D, solve2D_chunks, the
    stream (initial analysis and 2D pool), the particle sieve and `validate` all run at 64 unless given the key;
    the solve2D controller sets 1, not a removal, on stage 1, the seed pass and the terminal pass. The refusal of `topk_inpl > 1` with a refine mode without a table is removed:
    the key is read only where there is a table, and `solve2D` itself may run `refine=snhc_smpl`. Rulings 10 and
    section 3's off-state sentence are read with this: the byte identity of section 7 is that of `topk_inpl=1`,
    given explicitly where the default is 64.

## 5. Not in this plan

- **The exact angular marginalisation.** The posterior-weighted average over the whole circle is, per
  ring, a circular convolution of the particle's polar rings with the posterior, exact for every harmonic
  and continuous in angle, at a cost independent of `K`: a polar resampling of the particle's plane and its
  CTF-squared plane over the full lattice, an angular FFT per ring, the product with the posterior's
  harmonics, and a gridding of the polar samples into the Cartesian accumulators. It needs a polar geometry
  over every shell to the lattice's Nyquist (the calculator's covers the band) and a measurement of the
  two resamplings' blur against the direct insertion. It is the design to take if the `K`-point quadrature's
  cost at `K > m_max` or its surviving harmonics prove unacceptable on Hans's sets; it is not built here.
- The noise model: no change to `sigma2/`, to the residual at the committed pose or to the seed. Follow-ups
  that act on the loop's gain, in the order to try them once the insertion is in:
  1. Joe Caesar's early-stage shape: inflated sigmas at the rings beyond the FRC in stages 2 to 3, annealed
     to the measured curve (the postdocs' `sigma2_white` code is the instrument; the flat level everywhere
     is what collapsed their classes). Each ring's gain is the draw's factor times the ratio of the true noise
     to the assumed noise, so inflating the high rings early takes the top views' vote away from the rings
     with no angular signal.
  2. Posterior-weighted residuals (the sigma annealing proposal's section 7.5): the residual at the committed
     pose is biased low where the reference has power, which raises the effective temperature throughout.
  3. A reference-side Wiener factor per ring on the asymmetric angular harmonics of each class average from
     their cross-half correlation, which damps the residue at its source.
- The low-pass ladder (`lpstop` from the mask diameter, no data-driven cap below the FRC): a separate plan.
- `topk` and everything in `simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`,
  `simple_strategy3D_prob.f90`, `simple_matcher_3Drec.f90`; the stream (`src/main/stream`, the pool's
  `refine2D_distr`), the sieve, `make_cavgs`, the `cls_init=prev` seed pass, the carried sums
  (`cavg_state.bin`) and the contribution file; shift extras; class extras; a cone; the terminal pass.
- Larger-scale testing (ruling 4): fixture campaigns, real-data arms, floors, the real-data read; new
  fixtures, diagnostics or scripts. The unit assertions of Phase 1b are in.
- Skill edits (`simple-solve2d`, `simple-main-strategies-probabilistic`, `simple-main-class`): one line each
  is proposed in Phase 1 and applied only with approval.

## 6. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 baseline | none | Release and Debug builds, fast gate, the identity entries twice at one thread |
| 1 insertion | identical at `topk_inpl=1`, result-changing above | `test=identity` at one thread against Phase 0; fast gate; the existing unit suites; a smoke run at `topk_inpl=8` completes |
| 1b assertions | tests-only | each new check fails on its planted violation, then passes; `scripts/check_test_registry.py` |
| 2 handoff | none | the report and the runs named for Hans |

## 7. Acceptance

- **Identity** (the gate): `simple_test_exec test=identity` on `test=simulated_workflow` (`suite=1jxy`
  and `suite=6vxx`, each running solve2D and solve3D) at one thread with a fixed seed, base against the
  build at `topk_inpl=1`, exclusions declared from the two base runs of Phase 0 (an output that differs
  between the two base runs is excluded then, never later; the October run found the 6vxx entry
  irreproducible across thread counts, which one thread does not exercise). An item that cannot be made
  identical is restored and recorded; never rerun until it passes. This gate proves the off state only; it
  cannot detect a wrong sign, weight, index or accumulator on the enabled path (R1, point 3). That path is
  covered by the assertions below.
- **Assertions of the enabled path** (Phase 1b), in the existing suites, each planted with its violation
  first so that the check is shown to fail, then passing:
  - in `simple_eul_prob_tab2D_tester` (`suite=2d_probability_table`), on the quadrature selector:
    the weights of a particle sum to 1; coincident bins merge into one entry; the angle wraps at the grid's
    end; a row without a valid distance gives no entries (the committed pose alone); a flat row of `nrots`
    bins gives `K` entries at the keyed phase and, where `K` divides `nrots`, equally spaced with every harmonic
    up to `m_max` that is not a multiple of `K` within `1e-5` of the full-row sum's, otherwise (`K = 64` on 360
    rotations) the quantisation residual of the grid (12.5 % of harmonic 40; review R2, point 2) (and the offset's family draws
    the same phase on any thread); a row with all its mass in one bin collapses to one entry at the
    committed continuous angle at weight 1; a peaked row's response on `m <= m_max` is within the bound of
    ruling 15 of the full-row sum;
  - in `simple_pftc_inplane_tester` (`suite=continuous_in-plane`, which already exercises
    `gen_euclid_dists_bank`): the row's value at the committed rotation equals `gen_corr_for_rot_8_shifted`
    at the same pose for a non-zero particle-frame shift and an astigmatic CTF, within the kernel's
    single-precision spread measured on the pair;
  - in `simple_cavg_sums_tester` (`suite=class_average_carry_over`): one entry at weight 1 at the polished
    angle gives the same sums bit for bit as today's insertion; two entries at weight 0.5 at one angle equal
    one at weight 1 within rounding; the CTF-squared mass of a particle is the same for any set of weights
    summing to 1 (both accumulators take the weight); the sums with entries are bit-identical at 1 and 4
    threads.
- **The repository's standard gates**, unchanged: the fast gate (`scripts/run_fast_gate.sh`), the existing
  unit suites the touched folders name (`test=unit_pftc_align2D3D`, `test=unit_reconstruction
  suite=class_average_carry_over`, `test=polar_matching_bench suite=fills`), Debug build clean.
- **Smoke, not a test of the effect:** `solve2D` on the simulated fixture at `topk_inpl=8` in shared
  memory and `refine2D` with `nparts=2` at `topk_inpl=8` run to completion and print `PROB_TAB_INPL`
  with summed counts; the Debug build runs the shared-memory case once (bounds checks on the new arrays).
- **Reported, not gated:** the pass's wall time and `rt_cavg_interp_splat` at `topk_inpl` 1, 8 and 32 on
  the fixture (`bench=yes`), for Hans's cost reading; the montage of the fixture's averages at the last
  iteration of stage 6 and after the terminal pass, base and `topk_inpl=8` (paths in Progress).
- **The perturbation response, recommended by R1 point 4, run only if Hans asks:** one iteration of
  `refine2D` on the fixture from a set of references of which one is made rotationally symmetric and given
  a small asymmetric perturbation at one harmonic, at `topk_inpl=1` and at `topk_inpl > m_max`, and the
  amplitude of that harmonic in the resulting even and odd sums against its input. Attenuation below one
  supports the hypothesis of section 1; amplification refutes the temperature choice of ruling 14 for
  this purpose. A measurement for Hans's reading, not a gate.
- **Validation of the effect:** Hans's, after delivery, on the membrane-protein sets. The report names
  the command lines to run and the readout line to read: the mean `n_eff` and the collapsed fraction tell
  whether the top views' rows are flat (the quadrature spreads) or peaked (it collapses); `m_max` beside
  `K` tells whether the stage's harmonics lie off the multiples of `K`; the grid-snapping residual of review R2 remains.

## 8. Common rules and stop conditions

- No parameter but `topk_inpl` (ruling 2); registered in `simple_parameters.f90`, parsed, checked in the
  phases, in the UI for refine2D, solve2D and solve2D_chunks; `simple_args.f90` and
  `simple_ui_default_values.f90` are generated.
- Fix as you go, with a file table row before the edit; every file changed outside the table is listed in
  the report.
- Thread policy: per-thread work sized by the team in use (`omp_get_max_threads()` read before the
  allocation); no per-call allocation in the particle loop; byte identity judged at one thread only. The
  offset draw is keyed (ruling 15), so a particle's quadrature is the same on any thread, batch or part.
- Builds and runs through the agent kit (`~/agent_kit/agent_cycle.sh SIMPLE <action>`, `CLAUDE_CLI=/usr/bin/false`),
  never into the checkout's `build/`; long runs in the background.
- Compile-time policy: the new type's large components allocatable; lean `use ..., only:` imports; no new
  CMake targets or per-file options.
- Stops: the identity gate failing at `topk_inpl=1`; a file outside the table; a fixed ruling the code
  contradicts.
- Gates no ruling may waive: identity at `topk_inpl=1`; the fast gate and the existing unit suites.

## 9. Phases

- **Phase 0, baseline.** Release and Debug builds of `db20cde75` through the agent kit, the fast gate,
  machine and load recorded, the base build kept; `test=simulated_workflow` (both suites) twice at one
  thread with a fixed seed to establish the identity exclusions. Exit: the Progress entry holds the build
  records, the run directories and the exclusion list.
- **Phase 1, the insertion.** `prob_insert_inpl` and the quadrature selector (rulings 6, 15, 16) in
  `simple_eul_prob_tab_utils`, the family `RND_FAM_INPLQ` in `simple_rnd`; the per-thread bank in `s2D`;
  in `strategy2D_prob%srch` the row at the committed pose (ruling 7), the quadrature at the particle's own
  temperature (ruling 14) with the committed-bin substitution (ruling 16) and the fill of the particle's
  entries; `cavger_update_sums` consumes (ruling 8); the readout fields and the `PROB_TAB_INPL` line with
  `m_max`; the key through params, UI and the controller (ruling 10); the policy text (ruling 1) in
  `doc/policies/2D/solve2D_policy.md` (the stage controller section and the ownership section), the
  accumulation text of `doc/algorithms/refine2d_class_averaging.md` (its stale noise-model paragraph
  corrected while there, the sigma2 cancellation stated under the assumptions of section 2), and the
  sentence "the 2D route is untouched" of the refine3D policy pointed here. Exit: `test=identity` at
  `topk_inpl=1` passes against Phase 0; the fast gate, the unit suites and the Debug build pass; the smoke
  runs of section 7 complete and print `PROB_TAB_INPL`; the cost numbers and montage paths are in Progress.
- **Phase 1b, the assertions.** The checks of section 7 in the three testers, each written with its
  violation planted first (a wrong sign, a weight on one accumulator, an index off by one, a missing
  merge) and shown to fail, then the plant removed. Exit: every planted violation fails its check; the
  three suites pass; `scripts/check_test_registry.py` passes; the fast gate passes.
- **Phase 2, handoff.** The report `inpl_quadrature_insertion_2D_report.md` beside this note, in the
  order of `doc/how2s/how2write_agent_plans.md` section 8: what was learned, the files changed with counts,
  decisions taken under delegation, deviations, the identity evidence, and for Hans the command lines
  (`solve2D ... topk_inpl=8`, `topk_inpl=32`, and the `m_max` of his sets' stages) and the readout line
  to read. Exit: the report exists; the tree is left uncommitted with the changed files listed.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | `prob_insert_inpl`; the quadrature selector; readout fields and `PROB_TAB_INPL` | 1 |
| `src/utils/math/simple_rnd.f90` | `RND_FAM_INPLQ` | 1 |
| `src/main/strategies/search/probabilistic/simple_strategy2D_prob.f90` | the row at the committed pose, its `n_eff`, the quadrature fill, the committed-bin substitution | 1 |
| `src/main/strategies/search/simple_strategy2D_srch.f90`, `simple_strategy2D_alloc.f90` | the per-thread bank in `s2D`, its allocation and release | 1 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90` | the side channel per batch, handed to the sums; the readout counts; `m_max` | 1 |
| `src/main/class/simple_classaverager.f90`, `simple_classaverager_restore.f90` | `cavger_update_sums` takes the entries; inserts them in the particle's slot, each weight on both accumulators | 1 |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | `topk_inpl`, parse, range check | 1 |
| `src/main/ui/simple/simple_ui_refine2D.f90` | the key on refine2D, solve2D, solve2D_chunks | 1 |
| `src/main/solve/simple_solve2D_controller.f90` | the key on stages 2 to 6, deleted on stage 1 | 1 |
| `src/main/commanders/simple/simple_commanders_solve2D.f90` | the key deleted on the terminal pass; forwarded to the chunks' command line (`init_chunk_clustering`) | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D_tester.f90` | the quadrature assertions | 1b |
| `src/main/pftc/simple_polarft_mirr_tester.f90` | the row-at-the-committed-pose assertion (moved here from `simple_pftc_inplane_tester.f90`: its fixture has no CTF; this one has an angle-dependent CTF, per-ring noise weights and brute-force oracles) | 1b |
| `src/main/class/simple_classaverager_insert_tester.f90`, `simple_classaverager.f90`, `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90`, `doc/policies/test_environment_policy.md` | new: the weighted-entry and mass assertions on the private accumulation primitive, a submodule tester of the class averager, registered as suite `class-average in-plane entries` of `unit_reconstruction` (replaces `simple_cavg_sums_tester.f90`, which tests another type); the thread-count check of the sums is not unit-tested (section 7, Progress) | 1b |
| `doc/policies/2D/solve2D_policy.md` | the stage averages against the terminal and `make_cavgs` averages (ruling 1); the key | 1 |
| `doc/algorithms/refine2d_class_averaging.md` | the accumulation step; the stale noise-model paragraph; the cancellation under its assumptions | 1 |
| `doc/policies/3D/refine3D_policy.md` | one sentence: the 2D route's own insertion is here | 1 |
| `.github/skills/simple-solve2d/SKILL.md`, `simple-main-strategies-probabilistic/SKILL.md`, `simple-main-class/SKILL.md` | one line each, with approval only | 1 (if approved) |
| `doc/implementation_notes/planned/inpl_quadrature_insertion_2D_report.md` | the report | 2 |

## 10. Review R1 (2026-11-13) and what it changed

A review of the first version of this plan made five points and one implementation note; all were checked
against the code and hold.

1. A fixed midpoint-CDF quadrature of a flat row is an average over a cyclic group of order `K`, not over
   the circle: harmonics that are multiples of `K` survive at full amplitude, and the same `K` grid angles
   would be used for every flat-row particle. Changed: the keyed per-particle offset and the harmonic
   criterion `K > m_max` (ruling 15); the "no seed" claim withdrawn from section 1.
2. A grid quadrature discards the polished continuous angle, and a concentrated posterior would insert a
   well-determined particle at its bin centre. Changed: the committed-bin substitution and the collapsed
   case's fallback to today's insertion (ruling 16), now a contract line.
3. The identity gate at `topk_inpl=1` cannot show the enabled path is correct, and a smoke run that prints
   counts cannot either. The review asks for assertions in the existing suites: normalisation, merging of
   coincident bins, the invalid-row fallback and periodic wraparound; the angular-harmonic response of the
   `K`-point rule against the full-row sum; the narrow-posterior continuous-angle limit; non-zero
   particle-frame shifts and astigmatic CTFs; both accumulators, mass conservation and several threads.
   Changed: Hans restated ruling 4 (no larger-scale testing; unit assertions are in) and the list became
   Phase 1b.
4. "Unchanged loop gain" and "the physical posterior's expectation" overstated what a finite conditional
   quadrature gives, and the response of a posterior expectation to a small reference perturbation scales
   with the temperature, so the temperature change of ruling 14 changes it. Changed: section 1 states the
   reduced seed as a hypothesis and the gain as open; the perturbation measurement is in section 7 for
   Hans to order.
5. The sigma2 cancellation in the restored averages holds per half, per shell, for one iteration's sums,
   and not through the KB window's radial mixing, the merge of the halves before normalisation or the
   carry-over blend. Changed: the section 2 bullet and the algorithm doc's correction, under those
   assumptions.
6. The row is evaluated through the public `gen_euclid_dists_bank`, not the calculator's private helpers.
   Changed: ruling 7.

The review's closing points, that the benefit must survive the terminal hard pass and that "search
unchanged" holds for identical incoming references only, are in section 1 and the contract.

## 11. Review R2 (2026-11-13, after delivery) and what it changed

Two findings on the delivered code, both checked and both right.

1. `solve2D_chunks` accepted `topk_inpl` and dropped it: `init_chunk_clustering` builds the chunks' solve2D command
   line from scratch with an explicit key list, so every chunk ran at 1. The plan's premise (section 2) that a chunk
   copies the whole command line was wrong. Changed: `init_chunk_clustering` forwards the key when it is given;
   section 2, the file table and Progress corrected. Checked by running `solve2D_chunks` with the key (Progress).
2. The harmonic-cancellation statement was stronger than the quadrature gives: the points are snapped to the
   rotation grid, so a flat row's cancellation is exact only where `K` divides the number of rotations (for a phase
   that is not a tie). On 360 rotations `K = 64` retains 12.5 % of harmonic 40 at any phase, `K = 32` 13 % and
   `K = 16` 17 % of the harmonics up to 50, `K = 128` 6 % up to 80. `K > m_max` keeps the resolved harmonics off the
   multiples of `K` (no aliasing onto the surviving harmonics); it does not remove the residual of the grid.
   Changed: the statement is qualified in the policy, the algorithm note, section 1, ruling 15 and section 7 of
   this plan, the report and the selector's header comment; two assertions record the exact case (`K = 36`) and the
   residual (`K = 64`). The code is unchanged: the points stay on the grid, because a continuous angle inside a bin
   would assume the posterior uniform within the bin, which is wrong at a peak (ruling 16). A continuous
   interpolation for rows that are flat within a bin is the alternative if the residual matters on real sets.

## Progress

**Phase 0, baseline (2026-11-13, FWM-S116118, Apple M-series, gfortran 16.2, 24 cores shared with the maintainer's jobs).**
Base `db20cde75` built Release from a `git archive` into its own tree `inpl_base` (agent-run 1668): fast gate 15 of
15 passed. `test=simulated_workflow` for both suites on that tree at one thread (run 1669, ctest, registry command
lines): 6vxx passed in 792.5 s, 1jxy in 917.2 s; the run directories were copied aside (`inpl_ident/base_a`). A second
baseline run (agent-run 1683) and its comparison are in the Phase 1 entry. No exclusion list was needed (below).

**Phase 1, the insertion (2026-11-13).** All of section 9's Phase 1 items are in the working tree, uncommitted; files in
the report. Debug build of the working tree compiled on the first attempt and passed the fast gate 15 of 15 (run 1670),
the new suites included. Release build (run 1672): fast gate 15 of 15. Exit items:

- *Identity at `topk_inpl=1`.* New Release build, both suites at one thread (run 1678: 6vxx 759.5 s, 1jxy 883.5 s, both
  passed). `test=identity` between the run directories where the runs executed (`build_inpl_base/test_runs/...` and
  `build_release/test_runs/...`, runs 1681 and 1682 of the kit): 6vxx 177 files the same, 0 differing; 1jxy 214 files the
  same, 0 differing, no exclusion. (A first comparison of copies of the run directories reported the per-stage project
  files as differing: the tool strips each directory's own prefix from the records and the copies were elsewhere, so
  the paths of the original locations remained in the records; the comparison above is the one the tool is built for.)
- *The second baseline.* The base tree again (run 1683): 6vxx passed in 756.4 s, 1jxy in 879.5 s. `test=identity` between
  it and the new build's run directories: 6vxx 177 files the same, 0 differing; 1jxy 214 the same, 0 differing (runs
  1684, 1685). With the comparison above, the base is reproducible at one thread with the fixed seed and the exclusion
  list is empty.
- *After removing the upper bound of `topk_inpl` (Hans, 2026-11-13: `topk` has none),* the identity gate was run again on
  the final binaries: Release rebuilt, fast gate 15 of 15 (run 1686); `topk_inpl=0` is refused with its message and
  `topk_inpl=200` completes (runs 1687, 1688); both workflows passed (run 1689: 6vxx 778.7 s, 1jxy 907.0 s) and
  `test=identity` against the second baseline's run directories gives 177 and 214 files the same, 0 differing
  (runs 1690, 1691). The Debug tree of run 1675 predates the removal; the change is a deleted range check.
- *Review R2 (section 11).* Release rebuilt with the chunk forwarding and the grid-snapping assertions, fast gate 15 of
  15, `2D probability table` 126 checks (run 1692); `solve2D_chunks` with `topk_inpl=8` on the fixture (pruned first:
  runs 1693 and 1695) completed and the chunk's log holds 22 `PROB_TAB_INPL` lines (run 1696); Debug fast gate 15 of 15
  (run 1697); both workflows on the Release build passed (run 1698) and `test=identity` against the second baseline's
  directories gives 177 and 214 files the same, 0 differing (runs 1699, 1700), with no source file newer than run 1692.
- *The default of 64 (ruling 17; the report's section 8).* Release rebuilt, fast gate 15 of 15 (run 1705); `solve2D` without
  the key prints `topk_inpl 64` (run 1708) and `solve2D refine=snhc_smpl` completes (run 1709); the off state,
  `topk_inpl=1` against the base build at one thread with a fixed seed: 85 of 85 image and binary files identical,
  and 70 of 85 differ at 64 (a negative control); both workflows pass at the new default (run 1712).
- *Smoke runs at `topk_inpl=8`* on the extract-stage project of the 6vxx fixture (73 particles, 4 classes, `mskdiam=180`):
  shared memory, Release (run 1673) and Debug with bounds checks (run 1676, binary of run 1675), both completed with 22
  `PROB_TAB_INPL` lines; over two parts (`nparts=2`, run 1674) completed and the master printed the merged lines with 73
  particles (the parts' counts summed). Stage 1 prints none (the controller removes the key there); the control at
  `topk_inpl=1` (run 1677) prints none. Cost on this small fixture, matcher times summed over the 27 iterations of the six
  stages at four threads: total 2.22 s at 8 against 2.09 s at 1; class sums 0.32 s against 0.27 s; table fill 0.28 s
  against 0.34 s (the runs follow different trajectories). It cannot say more.
- *What the readout showed on the fixture,* for the maintainer's orientation: in the stage's first table iteration 93 %
  of the particles collapse to the polished pose at mean `n_eff` 1.0; later iterations 25 to 90 % collapse (mean
  `n_eff` 10 to 17 of 360 rotations, 1.4 to 2.3 entries per particle); `m_max` 50 at stage 2 to 71 at stage 6. At the
  physical temperature most rows of this simulated set are concentrated; the flat rows the plan is about would show as
  a large `n_eff` and a small collapsed share.

**Phase 1b, the assertions (2026-11-13).** Written in three existing testers and one new submodule tester (file table;
deviations in the report). Each was planted with its violation first: a build with six planted violations (agent-run
1671, Debug) failed exactly the two suites that hold the checks and no other: the lost CTF^2 weight failed four checks of
the class-sums suite; the reversed bank sign failed both frame checks (differences of 2.5 against a tolerance of 5e-5);
the missing merge, the ignored phase and the reversed weights failed ten checks of the 2D table suite. The plants were
reverted from backups (the three touched files equal their pre-plant content; two are equal to `HEAD`) and a clean
Debug build passed the fast gate 15 of 15 (run 1675): `class-average in-plane entries` 12 checks, `2D probability table`
122 checks (new tests included), `mirrored references` 4297 checks. Not unit-tested: the thread-count independence of
`cavger_update_sums` with entries (it needs the builder and a project; the bucket order is unchanged and a
particle's entries are applied inside its bucket in entry order). The key's forwarding through `solve2D_chunks`
was missing in the first delivery (review R2) and is added in `init_chunk_clustering`.

**Phase 2, handoff (2026-11-13).** The report is written; the tree is left uncommitted with the changed files listed
there.
