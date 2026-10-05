# Report: 3-D fractional-update balance refactoring

Date: 2026-10-05 (updated after Phase 3). Plan: `view_balanced_fractional_sampling_refactoring.md` beside this report (its Progress section
holds the detail). Run: autonomous run `balance_sampling` on dell, base commit 9ac9517d0, nothing committed. The
diffs of the three phases are in `/home/elmlundho/agent_runs/balance_sampling/review/phase_N.diff`; logs, scripts and
tables are under `/home/elmlundho/agent_runs/balance_sampling/scratch` (`logs/`, `cmp/`).

Terms: *fractional update* refines only a sample of the particles each iteration (`nsample` of them); a *sampling
unit* is one selected 2-D class; *FSC 0.143* is the resolution where the Fourier shell correlation of the two
independent half-map reconstructions falls to 0.143; *sweep* is the number of draws that visit every particle once.

## What changed

- **Phase 0 (`phase_0.diff`, plan only).** Baseline on a Release build of the base: `solve3D nsample=1000` with
  `balance=yes` (A, 3.84 and 3.93 Å final FSC 0.143) and with `partition=yes` (B, 3.98 and 3.84 Å), and
  `refine3D_states` on A1's final project (C, 10:59). The project has 50 selected classes, so `nclust` keeps its
  default of 20.
- **Phase 1 (`phase_1.diff`).** `balance=none|class|cavg` replaces `balance=yes|no` and the sampling use of
  `partition` and `clust_crit`. Sampling units are the selected 2-D classes; under `cavg` they are grouped by
  average-linkage clustering of the class averages on their correlation, and the quota is nested: equal over groups,
  equal over the classes of a group, lowest update count first. `clssmp.bin` carries the group as a fourth header
  value. Projection-direction sampling in `refine3D_states` and the `cls2D` cluster branch of `classify3D_refs` are
  deleted. The sweep comes from the unit table (`refine3D_states` initialization and iteration cap, `refine3D_auto`
  cap); a coverage table is printed before the first stage with warnings. `refine3D_states` holds one particle cohort
  per frequency block (internal `cohort_sampling`, `sample4update_rescore`). New unit tests cover the nested quota, the
  sweep, cohort rescoring, the file's group column and the trailing weight of a held cohort. Policies and the two
  sampling skills are updated.
- **Phase 2 (`phase_2.diff`).** Comparison runs; under the maintainer's ruling the per-state `updatecnt == 1`
  diagnostic added in Phase 1 was removed again, and the `class` equivalence was settled by a fixed-seed paired test.
  The plan moved to `completed/`.
- **Phase 3 (`phase_3.diff`, plan and report only).** Added by the maintainer: two more seeded case-A runs, one per
  binary, to test reproducibility under a fixed seed. Both same-binary pairs differ from iteration 1, so `solve3D` is
  not reproducible under a fixed seed and the base/new outcome difference is not evidence against the change.

## Decisions taken under delegation

All are recorded in the plan's Progress section (Phase 0, Phase 1, Phase 2 entries).

- Within a group, the remainder of the equal split goes to the classes whose particles have the lowest mean update
  count (ties: larger population, then index): deterministic for every distributed partition and rotating over
  iterations. The group level keeps the former equal-increment rule, so `class` allocates exactly as `balance=yes`.
- `solve3D` writes `clssmp.bin` for every fractional run (class units under `none`, because its initial greedy sample
  and the docked split checkpoint always drew from class units); the split assignment uses the run's `balance`, or
  `class` under `none`.
- `solve3D` with input volumes defaults to `class` and rejects `cavg` (random classes have no averages), as it
  rejected `partition=yes` before.
- `refine3D_states` and `classify3D_refs` leave particles inactive in `ptcl3D` out of their units.
- The `classify3D_refs` iteration cap was left unchanged (the plan names only `refine3D_states` and `refine3D_auto`).
- The `selection` program's `balance` was a UI entry its commander never read; its choices became `none|class` and
  `cavg` is rejected on its command line.
- Case C′ restores A1's run directory from `scratch/keep/` before it starts, because the kept project refers to its
  consensus map by absolute path.
- Phase 2 added four unseeded spread runs (two base, two `class`) as evidence when A′2 failed; they are not
  acceptance runs. Phase 3 copied the Phase 2 seeded-pair script under new run names (`seed_pair2.sh`) and compared
  the logs with a new script (`seed_logcmp.py`).

## Deviations from the plan

- One file-table row added in Phase 1 (`simple_oris_tester.f90`, named in the plan's tests but missing from the
  table) and one in Phase 2 (`simple_convergence.f90`, for the ruling).
- The per-state `updatecnt == 1` diagnostic of the plan's "Convergence" section was implemented in Phase 1 and removed
  in Phase 2 by the maintainer's ruling: under cohorts it always reads 0.
- The A′ acceptance as written (each A′ within the larger of the Phase 0 spread and 3 % of A's mean, i.e. ≤ 4.00 Å)
  failed for A′2 (4.29 Å); the phase stopped, and the maintainer's ruling replaced the criterion for A′ by a fixed-seed
  paired identity test of the sampling.

## Evidence against the baseline

| Case | Baseline (Phase 0) | New code (Phase 2) | Acceptance |
|---|---|---|---|
| A / A′ (`balance=yes` / `class`) | 3.84, 3.93 Å | 3.89, 4.29 Å | by ruling: sampling identical (below) |
| B / B′ (`partition=yes` / `cavg`) | 3.98, 3.84 Å | 3.93, 3.89 Å | pass (≤ 4.05 Å) |
| C / C′ (`refine3D_states`) | state 1 4.03 Å, state 2 10.20 Å, 10:59 | 4.08 Å, 12.09 Å, 21:50 | completes, report shown, cohort behaviour holds |

- Same-seed control (Phase 3): base build 3.89 and 4.47 Å, new build 4.41 and 3.93 Å, each pair differing from
  iteration 1 (`scratch/logs/p3_seed_logcmp.txt`).
- Paired fixed-seed test (`SIMPLE_SEED=20261005`, inherited by the workers): all 133 sampling rounds of all eight
  `solve3D` stages hold the identical particle set in the base (`balance=yes`) and new (`balance=class`) runs; final
  `sampled` and `updatecnt` agree for all 5,513 particles (`scratch/cmp/seed/compare_rounds.txt`). Final FSC 0.143:
  3.89 Å (base) and 4.41 Å (`class`); their alignments differ from iteration 1 with identical particles.
- C′: `% PARTICLES UPDATED SO FAR` per frequency block 27.374, 48.243, 62.369, 74.240, 83.998, 90.646, 94.444,
  97.270, 98.789, 99.359, 99.739, 100.000, constant inside every block; the coverage report gives a 9-block sweep
  against 12 planned blocks, visits per particle 1.41 to 12.0, no warning; every active particle updated. C′ runs 36
  iterations (4 × sweep) instead of C's 15, hence the doubled wall time; its projection-direction sampling drew 55.6 %
  of the particles per iteration whatever `nsample` was, C′ draws 27.4 %.
- Builds: Release and Debug show the base's linker warnings only (5 and 6), no compiler warnings; the fast gate passes
  (14/14) on the Phase 1 Debug build and the Phase 2 Release build; the unit areas of the touched code and
  `check_test_registry.py` pass.

## Open for the maintainer

- **De novo divergence (not caused by this change).** `solve3D` case A on beta-galactosidase sometimes goes wrong in
  its ab initio stages (stage-4 FSC 0.143 near 8–9 Å instead of 5.5–6 Å) and ends above the Phase 0 limit of 4.00 Å.
  All runs of this refactoring: base `balance=yes` 1 of 6 poor (3.84, 3.93, 3.89, 3.98, 3.89, 4.47 Å), new
  `balance=class` 3 of 6 poor (3.89, 4.29, 5.10, 3.89, 4.41, 3.93 Å), view-partition B and B′ 0 of 4, and 0 of 6
  identical runs in the earlier trailing_halfmap run. `class` samples the same particles as `balance=yes` in every
  round (Phase 2 paired test), and Phase 3 showed that `solve3D` is not reproducible under a fixed `SIMPLE_SEED` even
  on one binary: two same-seed runs of the base build, and two of the new build, part at iteration 1 in cFAR,
  orientation change and score, as the base/new pair did. The spread is in the multithreaded alignment numerics, and
  1 of 6 against 3 of 6 does not separate from chance. Open: the failure rate of `solve3D`'s random start on this data
  set, and whether a fixed seed should make it reproducible. A two-run spread is too small a baseline for acceptance
  on `solve3D` resolutions.
- **Coverage warnings.** With the `refine3D_states` iteration cap of 50 and three-iteration blocks, any sweep above
  about 17 blocks will trigger the "least-visited unit below one visit" warning; by design nothing is adjusted.
- **Pre-existing, not changed:** distributed non-probabilistic refinement draws its stochastic sample per partition
  with independent seeds, so the union is not one quota-exact sample; the `pose_policy=fixed` initialization phase may
  stop on the overlap target before a full sweep.
- **Old artifacts:** `clssmp.bin` files and `solve3D` run manifests recording `partition` from before this change do
  not read (release 4 rules); `.simple` projects are unaffected.
