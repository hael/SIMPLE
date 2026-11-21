# Probability-table scoring as matrix products: the report

The autonomous run `table_gemm`, 2026-11-21, on the Dell (Intel Core Ultra 9 285, 24 cores; gfortran 15.2.1, the
bundled OpenBLAS 0.3.34), at base `7a679f424`, carried out the plan
[prob_table_matrix_scoring.md](prob_table_matrix_scoring.md) in eight phases. The plan's Progress section holds the
detail and every path; the evidence is under `/home/elmlundho/agent_runs/table_gemm/scratch/keep/`, and the driver's
diffs per phase are `/home/elmlundho/agent_runs/table_gemm/review/phase_N.diff`.

## What was learned

In short: SIMPLE's probabilistic alignment now scores the dense 3D probability table (`refine=prob`, stages 3-5 of
`solve3D`) and the subspace representatives of the neighbourhood table's `state` mode (`refine=prob_neigh`, stages
6-8) as complex matrix products over blocks of particles, with the per-pair kernels' own distance conversion and
keyed in-plane draw. On beta-galactosidase (bgal: 4,967 particles, four parts of six threads) the dense table's score
loop costs 0.38-0.40 of what it did and the `state` mode's 0.53, at the base's resolution, and every fill is
identical at any thread count, batch size and part count. The neighbourhood's other candidates (the members of the
pooled cells, the `geom` cones) were measured on the new kernel, missed the bar of two thirds of the per-pair cost
by a wide margin, and stay per pair.

A *pair* is one particle scored against one reference over every in-plane rotation: two sums over the Fourier shells
per angular frequency bin, one inverse fast Fourier transform (FFT) over the rotations, the distances, and a draw of
the in-plane rotation from the pair's posterior at the particle's temperature, keyed on the particle and the
reference so that the same pair draws the same rotation on any thread. The *per-pair routine* is
`gen_posterior_val` (and its seed-bank form), the *block kernel* the new `score_block_posterior`, a *fill* the
parallel loop that writes a batch's table rows. *Thread-seconds* are the threads' busy time summed. FSC is the Fourier
shell correlation.

- **The bundled OpenBLAS's `cgemm` is not position-invariant in columns or shapes, but is in rows.** An element of a
  product gets other last bits when it sits in another column or the call has another shape (72 % of the elements of
  a 32 by 16 product differ from one-by-one products), but not when only its row moves within a full call (probes in
  `keep/blas_probe/`). The kernel therefore always runs one shape, 16 particles by 32 references, with absent
  particles as zero rows and each reference in the column of its place in the list. A pair's entry then depends on
  neither its block-mates nor the team, and every fills comparison of the run was identical
  (`polar_matching_bench suite=fills`).
- **The posterior draw amplifies last-bit differences.** The block kernel's losses agree with the per-pair routine's
  within 2.5e-7 to 7e-7 relative. That changes the best rotation of 2 of 320,000 pairs, the rate the plan had
  measured, but the rotation drawn at the particles' temperatures (about 420) in 70 of 320,000: the draw multiplies a
  loss difference by the temperature. The kernel's tester holds the draws at its sizes (0 differing); the
  neighbourhood benchmark's checks compare the best rotation and the losses, as the maintainer ruled.
- **The dense table gains what the benchmark promised where the sums are long, less where they are short.** Over the
  bgal run its score loop costs 0.38-0.40 of the per-pair routine's: 0.31 at stage 5 (2,500 directions, shells to 65),
  0.5-0.62 at stages 3-4 (1,000 directions, shells to 26-33). There the per-pair tail (distances and draw, unchanged)
  and the inverse transforms keep more of the cost, and blocks of 12 particles leave a quarter of the 16-row product
  empty. The dense fills take 1.4 times less wall time and the whole run 9 % less (12:57-13:18 against
  14:17-14:32).
- **bgal's `solve3D` reproduces itself bit for bit** at a fixed seed, across runs and across thread counts per part.
  The plan's map tolerance, "inside the base's two-run spread", was therefore zero, which no last-bits change can
  meet; the maintainer replaced it (FSC within one shell, correlation with the base's map at least 0.98). The adopted
  dense table's final map correlates 0.99897 with the base's at the base's FSC (3.75 Å at 0.143, 4.29 Å at 0.5); the
  base's own last-bits perturbation (three parts instead of four) correlates 0.9984.
- **The neighbourhood's per-particle candidates do not fit a fixed product.** At the matcher's batch a pooled cell is
  shared by 2.3-2.4 particles of a block (a cone's cell by 3.6-4.4), so the 16 by 32 product uses a tenth to a fifth
  of its pairs: the members cost 0.98-1.22 and the cones 1.20-2.66 times the per-pair routine at 24 threads, against
  0.32-0.48 for the representatives every particle scores. The per-pair routine is bound by memory on the members at
  24 threads (52 µs a pair against 14 at one thread: each particle reads its own references).
- **The `state` mode's adoption changes the cost, not the result, on the data tried.** On bgal its score loop costs
  0.53 of the base's and its fills 0.83 of the base's wall, and every statistic of every iteration and the final map
  are bit for bit those of the run before the adoption: no peak, pool or draw flipped. The 6vxx workflow profile
  writes the same metrics in all 20 executions as before the adoption, and the multi-state label-recovery route
  prints the same results as on the base build.
- **The derived reference memo costs memory.** The block kernel reads a bin-major copy of the reference memo, made
  for every reference: +1.7 GiB per part at 5,000 directions (peak 2.93 to 4.64 GiB) and +1.3 GiB per 6vxx
  execution, although the `state` mode reads only its 252 representatives of 5,000.
- **The benchmark's projection for the `state` mode was pessimistic.** The synthetic neighbourhood route projected
  0.62-0.68 of the per-pair score loop; bgal gave 0.53, probably because each of its four processes reads its own copy
  of the references, which makes the per-pair representatives dearer there (not measured).
- **A second production caller of the dense table was missed until the last phase.** The pose student's training
  (`learn_pose_student`) fills a dense table without asking for the kernel's memo and stopped at its first fill from
  Phase 2 on; the review of this report found it, a one-argument fix made it pass its test again (deviations below).

## What changed, phase by phase

| Phase | Kind | Files | What |
|---|---|---|---|
| 0, baseline | measurement | 1 (the plan) | base build kept; fast gate; fills identity on the base (48 comparisons); bgal twice (byte-identical) plus two perturbations; per-stage cost records; the 6vxx before-profile (14 of 20 pass, the base's known failures) |
| 1, the kernel | identical | 12 | the distance conversion and the draw made callable from a second submodule; the block kernel in the new `simple_polarft_block.f90` (175 lines) with its per-thread work and plan; the derived memo `memoize_ref_bins`; an explicit `cgemm` interface; the memory term; the kernel's tester in the fast gate; the benchmark sums the angular Nyquist bin and gains the `block_prod` layout. 6vxx identity against the base: 185 of 185 outputs identical |
| 2, the dense fill | last bits | 5 | `fill_tab_range` in blocks of one half, the block's rows scored at once; the matcher asks for the memo and counts it; two fills routes (seeded, no shift) |
| 3, the dense profile | last bits | 1 (the plan) | 6vxx after-profile: inside the spread, 17 of 20 executions identical to the base |
| 4, the neighbourhood measured | benchmark | 4 | `polar_matching_bench suite=neigh`: the `state` and `geom` candidate rules at bgal's late-stage sizes, the per-pair routine against three block paths |
| 5, the neighbourhood adopted | last bits | 4 | `state` mode's fill in blocks with the representatives on the block kernel; the memo through the pipeline constructor's option; fills routes `state` and `geom` |
| 6, the workflow profile | last bits | 1 (the plan) | 6vxx after-profile with both tables; the labels route; the summed wall-time guard applied to Phases 2 and 5 |
| 7, docs and report | docs, one fix | 6 | `refine3D_policy.md` section 8, the pointer in the layout benchmark's report, the README row, the plan, this report; the pose student's training asks for the memo |

Production code changed: `simple_polarft_{calc,core,corr,memo}.f90`, the new `simple_polarft_block.f90`,
`simple_linalg.f90`, `simple_matcher_memory.f90`, `simple_strategy3D_matcher.f90`, `simple_eul_prob_tab.f90`,
`simple_eul_prob_tab_neigh.f90`, `simple_strategy_pose_student.f90`. Test code: `simple_polarft_corr_tester.f90`,
the benchmark `simple_polar_matching_bench.f90`, its driver and its help text. Documents: the policy, the benchmark
report's pointer, the README row. In all 17 tracked files changed and three new files (the kernel's submodule, the
plan, this report).

## Decisions taken under delegation

Each is in the phase's Progress entry (the dead code in the phase journals and diffs), and in the plan's decision it
revises where it revises one.

- **The product's fixed shape** (Phase 1, revising the plan's layout and per-thread work): 16 by 32 whatever the
  block, the particle terms as the conjugated operand (the conjugation undone exactly in the transpose), the tile
  buffer one bin deep, one batched inverse transform of 32 per particle, the work and plan made with the calculator's
  other plans. Reason: the BLAS finding above; the plan's variable block and tile would have made the fills depend on
  the team.
- **The memo asked for by the caller** (Phases 2, 5 and 7, within the plan's memo design): a table's pipeline
  constructor builds the memo when its caller passes `block_memo`; the matcher, the pose student's training and the
  fills harness do. Reason: the probability-table tester builds a `state` table on a calculator without memoized
  references and never fills it; asking unconditionally failed the fast gate.
- **Four fills routes added** (Phases 2 and 5): the dense table without the polar shift map (seeded, with its bank)
  and without shifts, and the neighbourhood table's `state` and `geom` modes on a subspace map built as the builder
  builds it. Reason: the harness exercised neither the seeded dense path nor any subspace mode.
- **The neighbourhood benchmark's construction** (Phase 4, revising the route's construction): bgal's stage-8 sizes;
  synthetic references smooth over the directions, so that the peaks lie about a particle's previous direction as
  late in a refinement; the production candidate rules replicated and formed once from the per-pair entries; one
  kernel call per cell, every state's members in a list that depends on the cell alone (gathering several cells
  would make the entries depend on the block); the threads' busy time as the measure; the verdict on the sum of two
  repeats, written down before the runs; one more pass, reported only, of the per-pair routine over cell-ordered
  blocks.
- **The `state` mode's blocks ordered by the previous direction's cell** (Phase 5, revising the ordering): a block's
  particles then share pooled cells whose members one thread scores in turn (0.86-0.93 of the per-pair cost on the
  members in the benchmark).
- **Dead code removed**: the dense fill's per-pair reference scorers (Phase 2) and the neighbourhood's per-pair
  coarse pass (Phase 5) lost their last callers; `geom` keeps one seed-bank slot per thread.
- **References for the comparisons**: Phase 1's workflow identity used one 6vxx execution per build at the same seed
  and `test=identity`, because the identity tool compares two directories of the same execution; Phase 6 compared the
  labels route with the record of its 2026-11-05 profile in `multistate_posterior_state_draw.md`, because that
  profile's directory is on dell2, and ran it once on the base build for context.

## Deviations from the plan

- **Phase 2 stopped** at two exit items, the map against a zero spread and the dense bar on the median iteration
  (0.500 and 0.546 of the base against 0.5). The maintainer ruled the map tolerance (above) and restated the bar as
  the score loop summed over iterations and parts, since the median had weighed a short stage-3 iteration like a long
  stage-5 one. The same runs were judged again: 0.38 and 0.40, met. Later he restated the wall-time guard the same
  way (the fill wall summed, at most 1.05 times the base's), because the per-iteration guard is at the timing noise
  for 2-second fills: Phase 2's dense fills 0.672 and 0.720, Phase 5's `state` fills 0.832 and 0.835, both met.
- **The neighbourhood checks** compare the best rotation and the losses, not the drawn rotation (ruled after Phase 1
  flagged the draw's amplification).
- **A regression found in Phase 7.** Phase 2 made every dense fill require the memo and passed `block_memo` from the
  matcher, the only caller it searched; the pose student's training also fills a dense table and stopped with "stale
  or missing bin-major memo" from then on (reproduced on Phase 5's build). The fact-check of this report found it.
  The file was added to the plan's file table with Phase 7 and the call given `block_memo=.true.`; `test=pose_student`
  then passes every floor, as on the base build, with nine report rows a few thousandths apart (the merged map's FSC
  at 12 Å against the truth 0.74192 against 0.74189). No other caller fills a dense or `state` table without it.
- **Errors inside the parallel regions**: the block kernel stops (`THROW_HARD`) inside the fills' regions, as the
  per-pair kernels do, instead of raising a per-task status after the region as the plan's region rule asks.
- **Expected effects that did not occur**: the benchmark's 3.5 times per pair became 2.5-2.6 times over the dense
  score loop and 1.6-2.0 at short sums; two of the three neighbourhood paths missed their bar (the members at best
  equal to the per-pair routine); the 2026-10-21 report's concern that block order would end bit-reproducibility did
  not arise, because the draws are keyed.
- **Files outside the plan's first table**, each added to the table before its edit: the benchmark's driver and help
  text (Phase 1), the matcher (Phase 2), the benchmark again (Phase 5), the pose student's strategy and this report
  (Phase 7). Nothing was changed outside the table.
- **One memo for every reference in `state` mode**: the plan fixes the memo's form; it covers all references though
  the mode reads 5 % of them.

## Evidence against the baseline

| Gate | Result |
|---|---|
| Fast gate (15 entries, budget 60 s) | 15 of 15 after every edit; the kernel's tester 0.47-0.49 s |
| Kernel's tester against the per-pair routine | 0 drawn rotations differing (allowed 1 per case), drawn losses within 2.5e-7 relative (bound 1e-6), one thread against the gate's team of four identical in three repeats; 24 threads by hand identical in 30 repeats |
| Phase 1 identity, 6vxx at one thread | 185 of 185 outputs identical to the base build |
| Fills identity | 0 differing entries in every run: 48 comparisons on the base; 48, 72 and 72 in Phase 2; 96 per run in Phase 5 (before the adoption, then at 24, 24 and 4 threads) |
| Dense bar on bgal (summed score loop at most half) | 0.382 and 0.399 of the base (1,301 and 1,360 against 3,364 and 3,453 thread-s; with a third base run interleaved 0.377-0.399) |
| `state` bar on bgal (summed score loop at most two thirds) | 0.531 and 0.534 (309 and 311 against 591 and 574 thread-s) |
| Summed wall-time guard (at most 1.05) | dense 0.672 and 0.720; `state` 0.832 and 0.835 (for scale, the third base run: 1.025 dense, 0.991 `state`) |
| bgal map (FSC within one shell, correlation at least 0.98) | FSC 3.75 / 4.29 Å as the base (its two runs byte-identical, spread zero); correlation 0.99897 whole, 0.99886 in the mask; Phase 5's map identical to Phase 2's |
| Neighbourhood bar on the benchmark (at most two thirds of per-pair) | representatives 0.32-0.48 (passes); members 0.98-1.22 and cones 1.20-2.66 (negative results); checks inside the kernel's tolerances in all twelve cases; two runs agree |
| 6vxx profile (20 seeds, one thread) | after the dense table and after both: 14 of 20 pass as the base (the same six known failures); every gate row shifted 0 spreads but the truth correlation, +0.0002 (spread 0.029) |
| Labels route (three states) | every true state resolved, accuracy 0.999, adjusted Rand index 0.997, every floor met; the base build prints the same 41 result lines |
| Pose student (`test=pose_student`) | passes every floor after the Phase 7 fix, as the base build; stopped at its first fill before it |

The base build and the frozen builds compared are under `keep/` (`build_base`, `frozen/p1_kernel`, `p2_dense`,
`p4_neigh_v2`, `p5_state`, `p7_final`). Timings ran with nothing else heavy beside them, the load recorded with each.

## Left open for the maintainer

Reported, not gated:

- Memory: +1.7 GiB per part at 5,000 directions in `state` mode. A copy of the representatives alone would need about
  5 % of it; the plan fixed the memo's form, so this is a follow-up.
- Load balance of the fills falls with blocks as the unit of work: 0.99 to 0.90 (dense), 0.87 to 0.78 (`state`). The
  wall-time guard covers it.
- Whole-run wall: 9 % faster on bgal after the dense table; no further visible change after the `state` mode (its
  fills are about 50 s of a 13-minute run); 6vxx walls within the profile's scheduling noise.

Data- and machine-dependent:

- The identity of the fills rests on the BLAS keeping a row's bits wherever it sits in a full product, probed for the
  bundled OpenBLAS on x86. Another BLAS, or the Mac's, is not covered: run `polar_matching_bench suite=fills` there
  before relying on it. The kernel's tester in the fast gate checks one thread against four on every build.
- The dense gain depends on the shells summed (1.6-2.0 times at 26-33 shells, 3.2 at 65) and the `state` gain on how
  many processes read their own references.
- The high-level entries other than those run here (`pose_student`, the labels route, the 6vxx workflow) were not
  run; the nightly suite would cover them.

Follow-ups worth a plan:

- A representatives-only memo for the `state` mode (memory).
- The per-pair tail, the distance conversion and the draw (about 2.9 µs a pair at one thread), is now the larger part
  of a block's cost at short sums; the plan left it out of scope.
- The members of the pooled cells are memory-bound per pair at 24 threads; a kernel with a variable row count would
  need its own position probe before it could keep the fills' identity.
- The block kernel's errors inside the parallel regions (deviations above).

The plan and this report stay in `planned/` for the maintainer to move to `completed/` with the commit.

## Proposed lines for the pftc skill

For `.github/skills/simple-main-pftc/SKILL.md` (not applied; the maintainer edits skills). Under **Rules**, after the
memoisation order:

```
- The block kernel (`score_block_posterior`, submodule `_block`) reads a bin-major copy of the reference memo,
  `ft_ref_bins`: `memoize_ref_bins` after `memoize_refs`, outside any region; a generation counter stops the kernel
  on a stale copy. A table asks for it in its pipeline constructor (`block_memo`), and every caller that fills a
  dense or `state` table passes it (the matcher, `learn_pose_student`, the fills harness); a fill without it stops.
  It is released with the reference memos and at `kill`.
- The block kernel's product has one shape, `PFTC_BLOCK_NB` particles by `PFTC_BLOCK_NR` references (16 x 32), each
  reference in the column of its place in the list: the bundled OpenBLAS keeps a row's bits wherever it sits in a
  full call, not a column's or another shape's. Changing the shape, the column rule or the BLAS needs the kernel's
  tester and `polar_matching_bench suite=fills` again; without it the fills depend on the team.
- The per-pair and block kernels share the distance conversion and the draw (`euclid_dists_from_crvec`,
  `draw_posterior_rvec`, module procedures of the calculator's interface); a submodule never calls another
  submodule's private procedure.
```

Under **Where to start**, in the submodule list: "`_block` (the block kernel, its per-thread work and its batched
inverse-FFT plan)". Under **Tests**: "the block kernel's tester (`test_block_kernel`, polar correlation suite) holds
it to `gen_posterior_val`; `test=polar_matching_bench suite=block_prod` and `suite=neigh` time it against the per-pair
routine."

## The tree

`git status --porcelain` at the end of the run:

```
 M doc/implementation_notes/README.md
 M doc/implementation_notes/completed/polar_matching_layout_bench_report.md
 M doc/policies/3D/refine3D_policy.md
 M src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90
 M src/main/commanders/test/simple_polar_matching_bench.f90
 M src/main/pftc/simple_polarft_calc.f90
 M src/main/pftc/simple_polarft_core.f90
 M src/main/pftc/simple_polarft_corr.f90
 M src/main/pftc/simple_polarft_corr_tester.f90
 M src/main/pftc/simple_polarft_memo.f90
 M src/main/strategies/parallelization/simple_strategy_pose_student.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90
 M src/main/strategies/search/simple_matcher_memory.f90
 M src/main/strategies/search/simple_strategy3D_matcher.f90
 M src/main/ui/simple_test/simple_test_ui_highlevel.f90
 M src/utils/math/simple_linalg.f90
?? doc/implementation_notes/planned/prob_table_matrix_scoring.md
?? doc/implementation_notes/planned/prob_table_matrix_scoring_report.md
?? src/main/pftc/simple_polarft_block.f90
```

## Commit message

```
pftc: score the probability tables' dense rows and representatives as matrix products

The dense 3D probability table (refine=prob) and the neighbourhood table's state
mode (its subspace representatives) now score a block of up to 16 particles of
one half against 32 references per angular bin with one cgemm, transpose the
products, run one batched inverse FFT of 32 transforms per particle, and convert
and draw each pair with the per-pair kernels' own procedures and keyed variates
(new submodule simple_polarft_block.f90, score_block_posterior). The product has
a fixed 16 x 32 shape with each reference in the column of its place in the list,
because the bundled OpenBLAS changes an element's last bits with its column and
the call's shape but not with its row; the fills stay identical at any thread
count, batch and part count (polar_matching_bench suite=fills).

The kernel reads a bin-major copy of the reference memo (memoize_ref_bins),
asked for by the tables' pipeline constructors (block_memo, passed by the
matcher and by learn_pose_student) and counted by the matcher's memory estimate.
The per-pair kernels keep every other caller; the neighbourhood's pooled-cell
members and geom cones were measured on the kernel (polar_matching_bench
suite=neigh) and left per pair.

On beta-galactosidase (four parts of six threads) the dense score loop costs
0.38-0.40 of the per-pair routine's and the state mode's 0.53, the dense fills
1.4 times less wall time, the run 9 % less, with the base's FSC and a final map
that correlates 0.999 with the base's. Result-changing in the last bits; the
6vxx workflow profile is unchanged within its spread. Kernel tester in the
fast gate (unit_pftc_align2D3D).

Plan and report: prob_table_matrix_scoring.md and prob_table_matrix_scoring_report.md
(doc/implementation_notes).
```
