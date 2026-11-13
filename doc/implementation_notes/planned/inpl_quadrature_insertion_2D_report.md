# In-plane quadrature insertion in 2D: report

Written 2026-11-13 by the run of [the plan](inpl_quadrature_insertion_2D.md), base `db20cde75`. The first delivery is
committed as `d6cf4ae25`; section 8 (the default of 64) is uncommitted on top of it. Built and run on FWM-S116118 through the agent kit; run numbers are the kit's.

## 1. What was learned

1. **The gate holds: `topk_inpl=1` reproduces the base byte for byte.** Both simulated workflow entries (`1jxy`, `6vxx`)
   at one thread on the base Release build and on the new Release build, compared by `test=identity` over the
   registry's declared outputs in the directories the runs executed in: 6vxx 177 files the same, 1jxy 214 files
   the same, none differing, no exclusion declared (runs 1681, 1682). The second baseline run is identical to the new
   build on the same files (runs 1684, 1685), so the base is reproducible at one thread with a fixed seed and the
   exclusion list of Phase 0 is empty.
2. **The enabled path works end to end on the fixture and the checks that guard it can fail.** `solve2D` at
   `topk_inpl=8` completed in shared memory (Release 1673, Debug with bounds checks 1676) and over two parts (1674),
   printing the new readout line per iteration from stage 2 on; in the two-part run the master's line carries the
   parts' counts summed (73 particles). Six violations planted into the code were each caught by the check written
   for it (section 5); the clean tree passes the fast gate in Debug (1675) and Release (1672).
3. **The shifted-particle row is in the same frame as the committed shift, with an angle-dependent CTF.** The bank's row
   equals the direct sums of the particle shifted in its own frame and the shifted-reference kernel at the shift
   `R(theta_r) u` of each rotation `r`, within 5e-5, on a fixture with the astigmatic form of CTF and per-ring noise
   weights (`check_bank_particle_frame`). Reversing the bank's sign changes the difference to 2.5. This was the
   frame the whole insertion depends on.
4. **At the physical temperature most rows of the simulated fixture are concentrated.** In the first table iteration of
   stage 2, 93 % of the particles collapse to the polished pose at mean `n_eff` 1.0; later iterations collapse 25 to
   90 % of them at mean `n_eff` 10 to 17 of 360 rotations, 1.4 to 2.3 entries per particle. This is what the plan
   says the temperature does, and it is what the maintainer's sets will show for a symmetric top view: a large `n_eff`
   and a small collapsed share where the row is flat, a collapse where it is not. The 73-particle fixture cannot say
   more.
5. **The stage command line is one object reused across the stages.** Deleting the key at stage 1 (the stage without a
   table) would have removed it for every later stage. The controller therefore carries the user's value in the stage
   configuration and sets or deletes the key on every stage (found by reading, before any run).
6. **A copy of a run directory cannot be compared by `test=identity`.** The tool strips each directory's own prefix
   from the project records; the records hold the absolute paths of the place the run executed in. A first comparison
   of copies reported the per-stage project files as differing; the comparison of the original directories passes.

7. **`solve2D_chunks` dropped the key (found by review after delivery).** The UI accepted `topk_inpl`, but
   `init_chunk_clustering` builds the chunks' command line from scratch with an explicit key list and did not forward
   it, so every chunk ran at 1. The plan had assumed that a chunk copies the whole command line; it does not. Now
   forwarded (section 2); verified by a run of `solve2D_chunks` with the key (section 5).
8. **The harmonic-cancellation claim needed the grid qualification (review).** The points are snapped to the rotation
   grid. Exact equal spacing needs `K` to divide the number of rotations (for a phase that is not a tie); otherwise
   harmonic `m` retains up to about `m 2 pi / nrots`. Measured on 360 rotations: `K = 64` 12.5 % of harmonic 40 at any
   phase, `K = 32` 13 % and `K = 16` 17 % of the harmonics up to 50, `K = 128` 6 % up to 80. The policy, the algorithm
   note, the plan, the selector's comment and this report are qualified; two assertions record the exact case and the
   residual. The code is unchanged: moving the points off the grid assumes the posterior uniform within a bin, which
   is wrong at a peak (plan ruling 16); it is the alternative if the residual matters on real sets.

## 2. What changed, phase by phase

Phase 0: no source change; builds and baselines. Phase 1 and 1b: 23 files modified (876 insertions, 28 deletions) and two
files new (the tester below and the plan note); Phase 2: this report.

- **Key and chain** (`simple_parameters.f90`, `_parse.f90`, `_phases.f90`, `simple_ui_refine2D.f90`,
  `simple_solve2D_controller.f90`, `simple_commanders_solve2D.f90`): `topk_inpl`, default 64 (the parameter's, section 8), no upper bound as `topk` has none; on refine2D, solve2D,
  solve2D_chunks (forwarded to the chunks' command line in `init_chunk_clustering`); the stages with a table get the
  key's value, stage 1, the seed pass and the terminal pass get 1.
- **Selector, side channel, readout** (`simple_eul_prob_tab_utils.f90`, `simple_rnd.f90`, `simple_eul_prob_tab2D.f90`):
  `inpl_quadrature`; `prob_insert_inpl`; the keyed family `RND_FAM_INPLQ`; six fields of the readout record and the
  `PROB_TAB_INPL` line; the table carries the pass's counts to its readout.
- **The step** (`simple_strategy2D_prob.f90`, `simple_strategy2D_srch.f90`, `simple_strategy2D_alloc.f90`,
  `simple_strategy2D_matcher.f90`): the row at the committed pose through the shifted-particle bank and
  `gen_euclid_dists_bank`; the per-thread bank, allocated only above 1; the per-batch side channel; the counts.
- **The sums** (`simple_classaverager.f90`, `simple_classaverager_restore.f90`): optional entries on
  `cavger_update_sums`; a particle with entries is gridded once per entry at its weight on the numerator and the CTF^2
  sum, in entry order, inside its bucket.
- **Assertions** (section 4): `simple_eul_prob_tab2D_tester.f90`, `simple_polarft_mirr_tester.f90`, the new
  `simple_classaverager_insert_tester.f90` with its registration (`simple_commanders_test_class.f90`,
  `simple_test_ui_class.f90`, `test_environment_policy.md`).
- **Documents**: `solve2D_policy.md` (new section "In-plane insertion", including the difference in how the stage
  averages and the terminal and `make_cavgs` averages are generated), `refine2d_class_averaging.md` (the insertion; the
  stale noise-model paragraph corrected; the sigma2 cancellation stated under its assumptions), `refine3D_policy.md`
  (one sentence), the plan and this report.

## 3. Decisions taken under delegation

1. **Assertions of the enabled path live where their fixtures are.** The row check is in the mirror tester (it has an
   angle-dependent CTF, per-ring noise and brute-force oracles; the in-plane tester's fixture has no CTF). The weighted
   accumulation is tested on the private accumulation primitive, which needs a submodule tester of the class
   averager (`simple_cavg_sums_tester.f90` tests another type), registered as suite `class-average in-plane entries`
   of `unit_reconstruction`.
2. **A refine mode without a table reads no value of the key.** The first delivery refused `topk_inpl > 1` there; the
   refusal is removed (Hans's default of 64 for `solve2D` would have stopped `solve2D refine=snhc_smpl`).
3. **A quadrature in one rotation inserts the committed pose at weight 1 whichever rotation that is** (ruling 16 as
   written): a well-determined particle is summed as without the key even when the posterior's single rotation is
   not the drawn one.
4. **`m_max` is computed from the pass's own band limit** (`box_crop smpd_crop / kmax`) and the stage's `mskdiam`, and
   printed beside `topk_inpl`.
5. **The readout fields count rows evaluated;** a collapsed particle counts one entry and one collapsed row; a row
   without a valid distance counts nothing.
6. **The planted violations were applied to production files temporarily** (backups in `/tmp/inpl_plants`) and reverted;
   after the revert the three touched files equal their pre-plant content, two of them `HEAD`.

## 4. Deviations from the plan

- **Tester files differ from the plan's table:** the row assertion is in `simple_polarft_mirr_tester.f90`, not
  `simple_pftc_inplane_tester.f90`; the class-sum assertions are in a new `simple_classaverager_insert_tester.f90`, not
  `simple_cavg_sums_tester.f90`. Files outside the plan's table that were changed: `simple_eul_prob_tab2D.f90` (the counts
  into its readout), `simple_commanders_test_class.f90`, `simple_test_ui_class.f90`
  and `test_environment_policy.md` (the registration of the new suite).
- **The thread-count check of the sums is not unit-tested.** Plan section 7 asked that the sums with entries be
  bit-identical at one and four threads. `cavger_update_sums` needs the builder and a project; the bucket order is
  unchanged and a particle's entries are applied inside its bucket in entry order, but this was not verified.
- **The perturbation response of plan section 7 was not run** (it was to run only if the maintainer asked).
- **The fixture is the 73-particle 6vxx project of the baseline run,** not a larger set: the cost numbers below say
  nothing about real data.

## 5. Evidence against the baseline

| What | Result | Run |
| --- | --- | --- |
| Base build, Release, fast gate | 15 of 15 | 1668 |
| Base, both workflows, one thread | 6vxx 792.5 s, 1jxy 917.2 s, passed | 1669 |
| Working tree, Debug, fast gate | 15 of 15 | 1670 |
| Six violations planted, Debug | failed `unit_reconstruction` and `unit_pftc_align2D3D` only | 1671 |
| Working tree, Release, fast gate | 15 of 15 | 1672 |
| Smoke, shared memory, `topk_inpl=8`, Release | completed, 22 readout lines | 1673 |
| Smoke, `nparts=2`, `topk_inpl=8`, Release | completed, merged lines with 73 particles | 1674 |
| Working tree, Debug (clean), fast gate | 15 of 15 | 1675 |
| Smoke, shared memory, `topk_inpl=8`, Debug | completed, no runtime check failed | 1676 |
| Control, `topk_inpl=1`, Release | completed, no readout line | 1677 |
| New build, both workflows, one thread | 6vxx 759.5 s, 1jxy 883.5 s, passed | 1678 |
| Identity, base against new, 6vxx | 177 files the same, 0 differing | 1681 |
| Identity, base against new, 1jxy | 214 files the same, 0 differing | 1682 |
| Upper bound of `topk_inpl` removed (Hans); Release rebuilt, fast gate | 15 of 15 | 1686 |
| `topk_inpl=0` refused; `topk_inpl=200` runs, Release | exit 1 with the message; completed, 22 readout lines | 1687, 1688 |
| Final build, both workflows, one thread | 6vxx 778.7 s, 1jxy 907.0 s, passed | 1689 |
| Identity, second baseline against the final build | 6vxx 177 and 1jxy 214 files the same, 0 differing | 1690, 1691 |
| Review R2: chunk forwarding and the grid qualification; Release rebuilt, fast gate (`2D probability table` 126 checks) | 15 of 15 | 1692 |
| `solve2D_chunks` with `topk_inpl=8` on the fixture, pruned first (the first attempt stopped on micrographs without stacks, runs 1693, 1695) | completed; the chunk's log holds 22 readout lines | 1696 |
| Working tree, Debug, fast gate | 15 of 15 | 1697 |
| Final build, both workflows, one thread | 6vxx and 1jxy passed (898.7 s for the pair) | 1698 |
| Identity, second baseline against the final build | 6vxx 177 and 1jxy 214 files the same, 0 differing; no source file newer than run 1692 | 1699, 1700 |
| Second baseline run, both workflows, one thread | 6vxx 756.4 s, 1jxy 879.5 s, passed | 1683 |
| Identity, second baseline against new, 6vxx | 177 files the same, 0 differing | 1684 |
| Identity, second baseline against new, 1jxy | 214 files the same, 0 differing | 1685 |

The planted violations and what caught them (messages from the Debug run 1671):

| Plant | Caught by |
| --- | --- |
| the weight not applied to the CTF^2 sum | four checks of `class-average in-plane entries`: a quarter weight, two halves, weights summing to one, the CTF^2 mass |
| the bank shifts the particle the wrong way | both frame checks of `mirrored references` (difference 2.5, tolerance 5e-5) |
| coincident points not merged | the two-equal-rotations check and the one-rotation checks of `2D probability table` |
| the phase ignored | the flat-row check at v = 0.3 of `2D probability table` |
| collapsed rows not counted | the counts check of the side channel |
| the weights stored reversed | the entries check of the side channel |

Quadrature tolerances come from an emulation of the rule in Python on the same rows (flat K = 8: rotations
`14 + 45 (j - 1)`, harmonics 8 and 16 at 1.0; flat K = 32 at v = 0.7: harmonic 32 at 0.988, others up to 20 below
0.035; Gaussian row of ten rotations at K = 64: largest harmonic difference 0.0134 against a tolerance of 0.02;
masses 4/7, 2/7, 1/7 recovered exactly at K = 7).

Cost on the small fixture (73 particles, four threads, matcher times summed over 27 iterations of six stages,
`BENCH` files of runs 1673 and 1677): total 2.22 s at `topk_inpl=8` against 2.09 s at 1; class sums 0.32 s against
0.27 s; alignment search 0.50 s against 0.46 s; table fill 0.28 s against 0.34 s. The two runs follow different
trajectories; the numbers show only that the step is not expensive where most rows collapse.

## 6. Left open for the maintainer

**Reported, not gated.** The effect on the membrane-protein sets, which is the maintainer's testing:

```text
simple_exec prg=solve2D projfile=<project> mskdiam=<A> ncls=<n> topk_inpl=8 ...
simple_exec prg=solve2D projfile=<project> mskdiam=<A> ncls=<n> topk_inpl=32 ...
```

Read the line `PROB_TAB_INPL (topk_inpl K, m_max M)` after each `PROB_TAB2D_ASSIGN`: a large mean `n_eff` with a small
collapsed share means the rows of that iteration are flat and the quadrature spreads; a collapse means the posterior
at the physical temperature has one rotation, and the particle enters as without the key. `K` above `m_max`
keeps every resolved harmonic off the multiples of `K`, which survive at full amplitude for a flat row (the check of the
plan's review). The points are snapped to the rotation grid, so the cancellation of the other harmonics is exact only
where `K` divides the number of rotations; on 360 rotations `K = 64` retains 12.5 % of harmonic 40 and `K = 32` and
`K = 16` 13 % and 17 % of the harmonics up to 50 (finding 8). `K` above `m_max` is necessary, not sufficient.

**Open questions the runs do not answer:**

- whether the enabled update attenuates or amplifies a small asymmetric reference perturbation (plan section 7, not run);
- whether the benefit survives the terminal greedy pass, which sums every particle at its committed pose (the policy
  states this on purpose);
- the cost at `K` above 50 on a real set (each entry is one more gridding of the particle's plane);
- the default: left at 1.

**Follow-ups named in the plan, in the order to try them if the effect is small:** Joe Caesar's early-stage shape of
the noise spectrum beyond the FRC; posterior-weighted residuals; a reference-side filter on the asymmetric angular
harmonics. The exact angular marginalisation (a polar circular convolution per ring) is recorded in plan section 5.

**Skill lines, proposed and not applied** (plan section 5):

- `simple-solve2d`: "`topk_inpl > 1` sums a particle in the stages with a table at a K-point quadrature of its in-plane
  posterior (`strategy2D_prob`, `prob_insert_inpl`, entries of `cavger_update_sums`); the committed pose, the terminal pass
  and `make_cavgs` are unchanged (policy: In-plane insertion)."
- `simple-main-strategies-probabilistic`: "The 2D in-plane insertion is a step after the draw, not a table row: the row
  of the committed class at the polished shift (`memoize_shifted_ptcl`, `gen_euclid_dists_bank`) fills
  `prob_insert_inpl`; the readout record carries the `inpl_*` fields."
- `simple-main-class`: "`cavger_update_sums` takes optional entries (`ent_n`, `ent_e3`, `ent_w`): a particle with entries is
  gridded once per entry at its weight on both accumulators, else once at its committed angle."

## 7. The first delivery's commit message (committed as `d6cf4ae25`)

See the commit. Sections 1 to 6 describe that state, with the later review fixes included in it.

## 8. The default of 64 (Hans, 2026-11-13, after a run on beta-galactosidase)

A run of `solve2D` with `topk_inpl=64` on a beta-galactosidase set (5513 particles, 90 classes) gave class averages Hans
called almost perfect; the readout showed 1.1 to 1.8 entries per particle and 40 to 90 % collapsed, mean `n_eff` 1.0 to 1.4.
The default is now 64 (plan ruling 17). A first version set it in the solve2D commanders and pinned the stream, the
sieve and `validate` at 1; Hans then moved it into the parameter, so every 2D method inherits it:

- `parameters` declares `topk_inpl=64`; the interface shows `{64}` on refine2D, solve2D and solve2D_chunks; no commander
  or caller sets it;
- the solve2D controller sets `topk_inpl=1` (not a removal, which would now mean 64) on stage 1, the seed pass and the
  terminal pass;
- the refusal of `topk_inpl > 1` with a refine mode that has no table is removed, because `solve2D refine=snhc_smpl` would
  otherwise stop on its own default;
- the policy's defaults list and key paragraph, the plan (ruling 17) and this report say so.

Evidence for the first version (kit run numbers; binaries from run 1706, a Release build of `d6cf4ae25-dirty`, made by another session on the
tree that included these edits; the tree also holds that session's uncommitted 3D changes):

| What | Result | Run |
| --- | --- | --- |
| Release rebuild, fast gate | 15 of 15 | 1705 |
| `solve2D` without the key on the 6vxx fixture | prints `topk_inpl 64`, 22 readout lines, completed | 1708 |
| `solve2D refine=snhc_smpl` without the key | completed, no refusal, 22 readout lines | 1709 |
| Off state: new build at `topk_inpl=1` against the base build, one thread, `SIMPLE_SEED` fixed, 6vxx fixture | 85 of 85 image and binary files identical (class averages per iteration, even and odd, frcs, sigma2 state, carried sums) | direct runs, below |
| Negative control: the new build at the default 64 against the base | 70 of 85 differ | direct runs, below |
| Both simulated workflows at the new default, one thread | passed: 6vxx 636.8 s, 1jxy 925.2 s | 1712 |
| Static gates (`check_cli_chain`, `check_test_registry`, `check_descr`, `check_guide_imports`, `check_tester_linkage`) | pass | no run |

Notes on how the evidence was made:

- **The workflows are no longer byte-identical to the base at the default,** by design: they call `solve2D`, which now
  runs at 64, and their `solve3D` follows. Their gate is the registry's correctness gate, which passed. The byte identity
  of the off state is shown by the direct comparison above.
- **The direct comparison was run outside the kit's `run` action,** because the kit cannot pass an environment variable and
  a run needs `SIMPLE_SEED` fixed to be reproducible (the test runner fixes it for its own process; a first comparison
  through the kit's `run` without it reported every file differing, runs 1710 and 1711, and proves nothing). The two
  binaries were the base tree's `simple_exec` and a private copy of the new one (the other session rebuilds the shared tree),
  run with `SIMPLE_SEED=20260923`, `nthr=1`, in `inpl_smoke/offstate/{new,base,new64}`. Nothing was built or installed.
- **Not run:** `mini_stream_*`, `stream_chains`, `multistate_workflows`, `solve3D_addon`, `cont_refine3D_1jxy` and the other
  entries that call `solve2D`. With the default in the parameter the stream (initial analysis and 2D pool), the sieve
  and `validate` run at 64; their behavior there is untested.
