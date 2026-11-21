# Probability-table scoring as matrix products

Status, 2026-11-21: run `table_gemm` complete on the Dell (Phases 0-7), awaiting the maintainer's review: the
dense table and the neighbourhood table's `state` representatives score with the block kernel; the members and the
`geom` cones stay per pair. Report: `prob_table_matrix_scoring_report.md` beside this note. Planned and ruled
2026-11-20; rulings R1-R8 in Progress and the decisions.
Base: `7a679f424`. Supersedes: none.
Depends on: the bundled OpenBLAS staying on the calling thread inside an OpenMP team
(`SIMPLE_BLAS_TEAM_SAFE`, defined for every build since `4f7b5dcf0`; validated: yes, probed
2026-11-19 at up to 232 concurrent in-team callers); the measure-only benchmark `polar_matching_bench`
and its report (validated: yes, Mac 2026-10-21, Dell 2026-11-20, with the caveat of section 2 on its
objective).
Integration owner: `polarft_calc` (`src/main/pftc`: the derived memo and the kernel, in a new
submodule) and the two tables' fills, `fill_tab_range` (`simple_eul_prob_tab.f90`, the dense table)
and `fill_tab_subspace_range` (`simple_eul_prob_tab_neigh.f90`, the neighbourhood table); the layering
`pftc -> probabilistic table -> matcher` is unchanged. Queue: the continuous-view work
(`polar_continuous_view.md`) edits `simple_polarft_cont_calc.f90` and the matcher's view counters; this
plan edits neither. The multi-state validation of `multistate_sigma2_anneal_validation.md` runs on
dell2 and edits nothing here. Master is not frozen; the run does not sync with it.

## 1. Why

SIMPLE's probabilistic refinement scores every particle of a batch against a set of references, one
pair at a time (`gen_posterior_val`: the two Euclidean sums over the shells, one inverse FFT over the
rotations, the distances, the keyed in-plane draw). Two tables do it in a single-state `solve3D`:
the dense table (`refine=prob`) in stages 3-5, every reference; the neighbourhood table
(`refine=prob_neigh`) in stages 1-2 (`prob_neigh_mode=shc`) and 6-8 (`state`), and `refine3D_local`
(`geom`). On beta-galactosidase (`/data/bgal/tst`, 4,967 particles, four parts of six threads, box
256, point group D2):

- the dense table at 2,500 directions and shells to 65 (stage 5 of the run of 2026-11-15,
  `4_solve3D/BENCH/REFINE3D_BENCH_PART1.txt`, iteration 67): a table fill of 22.6 s per iteration and
  part, 79 thread-seconds of it in the score loop at 25.5 µs per pair and thread, 54 in the shift
  refinement, beside 27 s of partial reconstruction;
- the neighbourhood table in `state` mode at 5,000 directions and shells to 81 (stage 8 of the run of
  2026-11-06, `prev2/4_solve3D`, iteration 17): 309 references scored per particle, 252 of them the
  subspace representatives every particle scores (126·5,000/2,500, `NSPACE_SUB`), the rest the members
  of its three peak cells; a fill of 7.1 s, 18 thread-seconds of score loop at 42.8 µs per pair and 19
  of refinement.

The benchmark's `matrix_t` layout computes the two sums of every angular frequency bin as one
`cgemm` over a tile of references and a block of particles, transposes, and runs the inverse FFTs
batched. Measured against production's per-pair routine on 2,500 directions and two states (5,000
references, `pftsz` 182, shells 2-62; the Mac: 2026-10-21; the Dell, this machine: 2026-11-20,
`~/tmp/blas_review/pmb/polar_matching_bench.tsv`), µs per pair and thread:

| machine, threads | per_pair | matrix_t (best B, R) | speedup |
|---|---|---|---|
| Mac M3 Pro, 11 | 31.1 | 8.3 (32, 64) | 3.7 |
| Mac, 1 | 18.7 | 5.7 (32, 64) | 3.3 |
| Dell, 24 | 33.4 | 9.6 (16, 32) | 3.5 |
| Dell, 12 | 24.0 | 8.1 (16, 32) | 3.0 |
| Dell, 1 | 15.7 | 6.7 (32, 64) | 2.4 |

Smaller blocks cost little at the team: on the Dell at 24 threads B 4, 8, 16 and 32 give 11.3,
10.3, 9.6 and 10.1 µs at R 32; R 64 and 128 are no better. The sums shrink 3-4 times; the
distances and the draw (2.9 µs at one thread, the same code in every layout) stay. The dense table's
references are shared by every particle of a batch, and so are the neighbourhood table's subspace
representatives (about four fifths of its pairs in `state` mode): both are matrix products as they
stand. The members of each particle's peak cells, and `geom`'s cones, differ per particle: blocks of
particles that share cells must be formed first, and whether that pays is measured before it is
adopted (Phase 4).

## 2. Contract

- A row entry keeps its meaning: for every (particle, reference) the Euclidean distances over all
  in-plane rotations from the sums `sum_k w_k [FT(CTF2) FT(REF2) - 2 FT(X.CTF) FT(REF)*]` over the
  angular bins 1..pftsz+1 (the Nyquist bin included, as production has summed it since `9933db06c`,
  2026-11-03), the particle's norm and temperature, and the keyed draw
  `u(scope, particle, INPL, reference)` (`gen_posterior_val`). The block kernel changes only the order
  of the arithmetic of the sums. The distance conversion and the draw are the same compiled procedures
  for both kernels (decision 14). The draws are keyed on the particle and the reference since
  `prob_align_batch_pipeline`, so the block order changes no draw given equal distances.
- The neighbourhood table scores the same references per particle as today, in the same order of
  record: the representatives first, in subspace order, then the members of the pooled cells; the
  peak selection and the pooling read the same distances in the same order. Only the arithmetic of the
  sums changes.
- Precision: single, as the per-pair kernel already sums (`sigma2_wk`, `euclid_shell_sum`).
- What changes: the last bits of the distances. Measured on both machines, `matrix_t` against the
  benchmark's replica of the per-pair sums: at most 4.0e-7 relative, the best rotation differs in 2 of
  320,000 pairs. The benchmark's replica and matrix layouts still drop the angular Nyquist bin, which
  production sums since 2026-11-03: that is why on the Dell every benchmark vector differs from
  production's by up to 2.7e-3 absolute. It is a formula difference in the harness, not rounding.
  Phase 1 brings the harness to production's objective; the gate on the new kernel is its tester
  against production's own per-pair routine (decision 11), never the harness's replica. Each
  adoption is result-changing in the last bits, judged by profile, never by byte identity against the
  base.
- Thread identity holds as today: the new fills at one thread and at any team give the same rows,
  batches and parts included (`polar_matching_bench suite=fills`).
- Off-state: there is no switch. The per-pair routines stay for their other callers (the shift
  refinement, the `shc`/`snhc` neighbourhood search, the 2D tables, the cross-half diagnostics, the cc
  objective) and for any neighbourhood path Phase 4 does not adopt; their arithmetic does not change,
  and Phase 1 keeps every production output byte-identical.
- No new parameter. No environment variable: the BLAS inside the team is single-threaded by the
  build's contract.

## 3. Decisions

Ruled 2026-11-20 by Hans: decisions 1 (re-ruled the same day on corrected facts: the first menu said
only the dense table had a cost record), 3, 4, 6, 10, 11, 12, 13 by explicit choice, 16 and 17 by the
re-ruling of decision 1 ("dense, then neighbourhood, adopted only past a bar fixed before the
measurement"); the others follow from them or are fixed by measurement. Never renumbered.

1. Scope: first the dense 3D table's all-reference score loop in `fill_tab_range`, in its three paths
   (`process_particle_with_shift`, `process_particle_with_map`, `process_particle_without_shift`);
   then the neighbourhood table as decisions 16-20 say. `objfun=euclid` (both tables require it). The
   2D tables and the `shc`/`snhc` neighbourhood search stay on the per-pair routines.
2. Layout: the benchmark's `matrix_t`: per block of B particles of one half and tile of R
   references, one `cgemm('C','N')` per angular bin 1..pftsz+1, depth 2·nk (both sums fused), the
   products transposed so that each pair's inverse FFT reads contiguous memory, then R·B transforms
   in one batched c2r. A tile of consecutive references is read in place from the bin-major memo; a
   tile given as a list of references (the representatives, a cell's members) is gathered per bin
   into the thread's tile buffer.
   Revised in Phase 1 (rulings R2): the product has a fixed shape, 16 particles by 32 references,
   whatever B and the tile (an absent particle is a zero row, a partial tile's missing references are
   zero columns), with the particle terms as the conjugated operand and each reference in the column
   of its place in the list (tiles aligned to the list's start). The products are the conjugated sums;
   the transpose conjugates them back, exactly. The inverse FFTs run as one batched c2r of the 32
   transforms of a particle's tile, once per particle of the block. A full tile of consecutive
   references is read in place; any other tile is gathered bin by bin. Reason: the bundled OpenBLAS
   gives an element of `cgemm` other last bits in another column or in a call of another shape
   (probed on the Dell: up to 72 % of the elements), but not in another row of a full call; decision
   3's B depends on the team and the batch, so in the benchmark's layout a pair's last bits would
   depend on its block-mates and the fills could not be identical across threads, batches and parts
   (section 2).
3. Tiles: R = 32. B = min(16, max(4, nb_half / (2·nthr))), nb_half the batch's particles of the
   half, so that a batch gives the team at least two blocks per thread. Origin: the Dell sweep of
   section 1 (B 8 keeps 3.2 of the 3.5 at the team; B 4 is the floor, 3.0); frozen.
4. The unit of work is a block. Each fill keeps one parallel region per batch, dynamic over blocks;
   the batch's particles are ordered by half and cut into blocks that never mix halves (as the
   benchmark's `make_blocks`). A task processes its block whole: each particle's seed, the block's
   scores into the particles' rows or candidate buffers, then each particle's remaining steps in the
   order its path uses today. Per-thread state stays per thread; the dense fill's comment "a particle
   is processed whole" becomes "a block is processed whole". The dense fill's per-thread copies of a
   row (`dists_refs`, `inpl_refs`) are replaced by reads of `rows_batch(:,i)`.
5. The seed bank (the particle shifted once by its seed, `memoize_shifted_ptcl`) holds the block's
   particles per thread, `seed_bank(:,:,B,nthr)`; a particle without a seed shift uses its memo
   (`ft_ptcl_ctf`, `ft_ctf2`). The block kernel takes, per particle, which of the two it reads.
6. The bin-major memo `(2·nk, nrefs, pftsz+1, parity)`, `conj(FT(REF2))` and `FT(REF)` per bin, is
   derived from the finished `ft_ref`/`ft_ref2` by one layout transform (parallel over references,
   outside any region, as every memo), only when a table that uses the kernel asks for it: its
   `new_pipeline` calls it after the matcher's `memoize_refs`. It has no logic of its own: the mirror
   rule of `memoize_refs` is already in its source. It is released with the other reference memos and
   at `kill`. Memory: one more term of 2·nk·(pftsz+1)·nrefs·2·8 bytes (about 355 KB per reference at
   the benchmark's sizes, 0.9 GB at 2,500 directions and one state), counted by `matcher_memory_terms`
   only when such a table runs.
7. Per-thread work, in the calculator's existing per-thread pattern: the particle terms
   `(2·nk, B, pftsz+1)`, the tile `(2·nk, R, pftsz+1)` for gathered tiles, the products
   `(R, B, pftsz+1)`, the transposed products `(pftsz+1, R·B)` and the FFT output `(nrots, R·B)`,
   sized at `pftc%new` by the team and the largest B and R; a partial block or tile runs the full-size
   plan on zero-filled rows. One batched c2r plan of R·B transforms per calculator, made and destroyed
   with its other plans.
   Revised in Phase 1 (rulings R2, with decision 2's revision): the tile buffer holds one bin,
   `(2·nk, R)`, gathered just before its product; the products are `(16, 32, pftsz+1)`; the plan
   has R transforms and runs once per particle of the block, so an absent particle costs no
   transform. The work and the plan are made in `allocate_memo_workspace` (called by `new`) with the
   calculator's other plans, under the same planning flags and wisdom, and released in
   `kill_memo_workspace`.
8. BLAS: an explicit `cgemm` interface in `simple_linalg` beside `sgemm`'s; the call inside the team
   stays on the calling thread by the build's contract (both OpenBLAS flavours).
9. Acceptance class: result-changing in the last bits (section 2). Section 6 holds the numbers.
10. Keep rule for the dense table, frozen before the run that judges it: on the bgal run, over the
    dense stages' iterations, the median per-iteration `match3D fill score loop (thread-s)` of the parts
    at most half the base's, and no iteration's `match3D fill region wall time (s)` above the base's by
    more than 5 % plus the base's two-run spread. The wall-time guard also covers the load balance of
    blocks as the unit of work. Origin: the benchmark's 3.0-3.5 at the team less the draw's fixed
    share (guessed from the benchmark's decomposition; the draw is about a third of the kernel's cost
    at 24 threads).
    Restated by rulings R6 (2026-11-21, applied first in Phase 2): per adopted run, the dense score
    loop summed over every dense iteration and part at most half of the base's sum (the mean of the
    base runs); the wall-time guard unchanged. The median per iteration weighed a stage-3 iteration
    (about 10 thread-seconds) like a stage-5 one (about 70); the bar itself stands.
    The wall-time guard restated by rulings R8 (2026-11-21, applied first in Phase 6): per adopted run,
    the table's `match3D fill region wall time (s)` summed over its iterations and parts at most 1.05
    times the base's sum (the mean of the base runs). The base `solve3D` reproduces itself bit for bit,
    so its two-run spread is about zero and the per-iteration guard was 5 % of a 2-second fill for the
    short `state` iterations, at the timing noise: a third base run fails it once (2.24 s against
    2.18 s). Phase 2's verdict stands in both forms (Phase 6's entry).
11. The kernel's tester runs in the fast gate, in the polar correlation suite
    (`run_all_polarft_corr_tests`, `unit_pftc_align2D3D`), so every build holds the two
    implementations of the score together. It compares the block kernel with `gen_posterior_val` on
    a synthetic calculator of 200 directions at one and two states and 32 particles, with and without
    a seed bank, with consecutive and with listed tiles: the drawn rotation equal in all pairs but at
    most max(1, 2e-5 of the pairs) (measured rate 6e-6: 2 of 320,000, on both machines), the
    distances within 1e-6 relative (measured 4.0e-7 on both machines, times 2.5), and the block kernel
    bit-identical between one thread and the gate's team of four in three repeats (measured: identical
    in five repeats at 24 threads). Its cost is held under one second at one thread (the gate's 60 s
    budget).
12. Phase order: the kernel and the memo first (no production caller, outputs unchanged), the dense
    fill second, its profile third, then the neighbourhood (decisions 16-20). If Phase 0's bgal record
    shows the dense score loop below a quarter of the dense fill in every dense iteration, the run stops
    before any edit and reports the profile.
13. The kernel, its per-thread buffers and its plan live in a new submodule of `polarft_calc`,
    `simple_polarft_block.f90`; the type gains the components and the interfaces only. The
    correlation submodule does not grow beyond decision 14.
14. One distance conversion and one draw for both kernels. `euclid_dists_from_crvec` and the body of
    `draw_posterior_crvec` are private procedures of the correlation submodule today; a second
    submodule cannot call them (gfortran before 16 does not link it; `check_submodule_linkage.py`
    fails the gate). Phase 1 makes the conversion a separate module procedure declared in the
    calculator's interface, and splits the draw into one over an explicit real vector, with
    `draw_posterior_crvec` its one-line caller on `crvec1(ithr)%r`. No arithmetic changes: Phase 1's
    byte-identity gate covers it.
15. The benchmark's replica and matrix layouts sum the angular Nyquist bin as production does
    (Phase 1), so their timings and checks describe today's objective. Reported, not gated: the
    replica's difference from production falls from 2.7e-3 absolute to the rounding level.
16. The neighbourhood's scope, in three paths, each measured in Phase 4 and adopted in Phase 5 only
    past decision 17's bar: (a) `state` mode's coarse pass, every particle of a block against the
    subspace representatives of every active state (`find_peak_subspaces`), a listed tile of the same
    references for every particle; (b) `state` mode's member pass (`evaluate_neighborhood`), the cells
    of the block's pooled neighbourhoods, each cell's members against the sub-block of particles that
    pooled it; (c) `geom` mode's cones (`evaluate_cone`, `refine3D_local`), the reachable cells of the
    block's cones, each cell's members against the sub-block that reaches it, with today's cone test
    deciding which entries are recorded (the pairs outside a particle's cone are wasted work and count
    against the path). The `shc`/`snhc` search (a hill climb that stops at the first improvement per
    particle) is out of scope.
17. The neighbourhood's bar, frozen now, before Phase 4's measurement, per path: adopted only if (i)
    in the benchmark's neighbourhood route at the late-stage sizes (5,000 directions, shells 2-81,
    252 representatives per state, 3 peaks; one and two states; 24 threads) the path's score
    thread-seconds over the pairs a particle needs, wasted pairs included, are at most two thirds of
    the per-pair routine's (1.5 times or better), and (ii) after adoption, on the bgal run's `state`
    iterations, the median `match3D fill score loop (thread-s)` of the parts is at most two thirds of
    the base's, with decision 10's wall-time guard. Origin: guessed; below 1.5 times a second scoring
    path in the neighbourhood table is not worth its upkeep.
    Part (ii) restated by rulings R6 (2026-11-21): the `state` iterations' score loop summed over
    iterations and parts at most two thirds of the base's sum (the mean of the base runs). Its
    wall-time guard is decision 10's as rulings R8 restate it: the `state` iterations' fill region wall
    summed over iterations and parts at most 1.05 times the base's sum.
18. Blocks in the neighbourhood fill follow decision 4. For paths (b) and (c) the batch's particles
    are ordered by half and then by the subspace cell of their previous direction, so that a block's
    pooled cells and cones overlap; path (a) needs no ordering. The ordering changes which particles
    share a block, never what a particle scores.
    Revised in Phase 5 (rulings R2, the ordering key): the `state` mode's blocks, which adopt path (a)
    alone, are ordered by half and then by the cell of the previous direction too. The members of the
    pooled cells stay on the per-pair routine, and in this order one thread scores the members of a
    block's particles, which share cells, in turn: in Phase 4's route the per-pair routine costs 0.86-0.93
    of its particle-order cost on the members at 24 threads (measured before this phase's bgal runs).
19. The benchmark gains a neighbourhood route (Phase 4): the production `state` and `geom` candidate
    rules on the late-stage geometry (5,000 directions, D2 and C1), the per-pair routine and paths (a),
    (b), (c), with their pair counts, wasted pairs and thread-seconds, and the check of each path
    against the per-pair routine on the same pairs (decision 11's tolerances).
    Restated by rulings R7 (2026-11-21, applied first in Phase 4): the check compares the best rotation
    and the distances with decision 11's tolerances, the draw made the best rotation by a temperature
    scale of 1e25 as in Phase 1's sensitivity check; at the particles' temperature the draw amplifies
    last-bit differences, so the drawn rotations' differences are reported, not gated. Decision 11's
    fast-gate tester is unchanged.
    Built in Phase 4 (rulings R2, the route's construction): `suite=neigh` of `polar_matching_bench`
    (`polar_neigh_bench`) at bgal's stage-8 sizes (box 256 at 1.275 Å cropped to 162, mask 180 Å, lp
    4.03 Å, shells 2-81, `pftsz` 140), solve3D's 252 representatives, the production defaults of
    `npeaks` (3) and `athres_cont` (10°), and the matcher's batch at 24 threads (1,200 particles, blocks
    of 12 by decision 3). The references are smooth over the directions (random anchors at the subspace
    spiral, blended by the point group's view distance) and each particle is its previous direction's
    reference, rotated in plane, under its CTF and noise, so that the peaks lie about that direction as
    late in a refinement. The candidates are formed once by the replicated production rules from the
    per-pair routine's entries. Path (b) and (c) score one call per cell, every state's members in a
    list that depends on the cell alone: each reference keeps its column (decision 2's revision), so the
    entries do not depend on the block; gathering several cells into one list would make them depend on
    it. The timing is the threads' summed time in the scoring region, production's `score loop
    (thread-s)`.
20. A path that fails decision 17 (i) is recorded as a negative result with its numbers and left on
    the per-pair routine; Phase 5 adopts only the paths that pass. If none passes, Phase 5 records that
    and edits nothing, and Phase 6 is skipped with a Progress entry.

Defaults taken without a question: the memo component `ft_ref_bins`, the kernel
`score_block_posterior`, the submodule name of decision 13, the tester's sizes of decision 11, decision
17's 1.5 times, decision 18's ordering key, the run name `table_gemm`.

## 4. Not in this plan

The 2D tables (`simple_eul_prob_tab2D.f90`), the `shc`/`snhc` neighbourhood search, the shift
refinement and the seed estimation (`refine_best_refs`, `refine_best_neighbors`, `estimate_shift_seed`,
the search objects), the cc objective (`gen_corrs`), the posterior diagnostics' cross-half scoring,
the arithmetic of the distances and the draw (their 2.9 µs per pair is the next ceiling and a
separate note; decision 14 only makes them callable from the new submodule), any change to
`ft_ref`/`ft_ref2`, the continuous-view calculator, the matchers, the UI and the parameters, the
skills (the report proposes the pftc skill's lines; the maintainer applies skill edits). Out of
bounds: every file not in the table.

## 5. Change types and gates

| Phase | Kind | Gate |
|---|---|---|
| 0 | measurement | base build kept; the numbers of section 6 recorded before any edit |
| 1 | identical (the shared tail moved, the kernel and the derived memo without a production caller, the harness's objective) | fast gate with the new tester (decision 11); `unit_pftc_align2D3D`; `test=identity` of `simulated_workflow_6vxx` at one thread against the base build |
| 2 | result-changing, last bits | `suite=fills` (thread identity of the new dense fill); bgal twice: decision 10, map and FSC within the tolerance of section 6 (rulings R5) |
| 3 | result-changing, last bits | `simulated_workflow_6vxx` 20-seed profile after against before: inside the spread, no gate row failing |
| 4 | measurement (benchmark only) | decision 17 (i) per path; the benchmark's checks inside decision 11's tolerances (the best rotation and the distances, rulings R7); fast gate |
| 5 | result-changing, last bits | `suite=fills`; bgal twice: decision 17 (ii), map and FSC within the tolerance of section 6 (rulings R5) |
| 6 | result-changing, last bits | the 6vxx profile after both against before; `multistate_workflows suite=labels` once, reported against its profile of 2026-11-05 |
| 7 | docs | fast gate |

## 6. Acceptance

- Kernel tolerances: decision 11 (measured, both machines; frozen).
- Speed: decision 10 for the dense table, decision 17 for the neighbourhood (guessed; frozen now).
- Thread identity: `suite=fills` reports zero differing entries at one thread against the team, in
  batches and in three parts, for every route (measured: the base passes it; recorded in Phase 0).
- bgal maps (rulings R5, 2026-11-21, applied first in Phase 2): the adopted run's final FSC at 0.143
  and at 0.5 within one Fourier shell of the base's, and the correlation of its final map with the
  base's final map at least 0.98, over the whole map and within the mask (as Phases 0 and 2 measured
  them). Origin: the ruling. The base `solve3D` reproduces itself bit for bit under a fixed seed, so
  the two-run spread the plan first named is zero and cannot judge a last-bits change; for scale, the
  base's own last-bits perturbation (three parts of eight threads instead of four of six, Phase 0)
  correlates 0.9984.
- 6vxx: the entry's gate rows inside the before-profile's spread (`s` of the profile tool), 20 seeds
  at one thread; a failure rate not above the before-profile's.
- Reported, not gated: the per-iteration fill records of the bgal runs (the `match3D fill ...` lines
  of `BENCH/REFINE3D_BENCH_PART*.txt`) by stage and table, before and after each adoption; the fills'
  load balance; the memory estimate's new term; the harness's replica against production after
  decision 15; the bgal run's whole wall time before and after.
- The base's known failures are recorded in Phase 0 and left (the nightly triage note lists them).

## 7. Common rules and stop conditions

No new parameter (ruled). No floor or tolerance relaxed; a failed tolerance is a diagnosis. Fix as
you go, with a file-table row before the edit. Every parallel region added or changed follows the
rules of the Fortran guide (per-thread arrays sized by the team at the region, no strings inside,
a per-task status raised after) and of the pftc guide (new scratch in the calculator's per-thread
buffers, never a per-call allocation in the score loop). One thread for every identity comparison.
Builds for before/after timings from one snapshot of the tree, interleaved, two repeats each way, with
the load recorded. Rulings at launch, applied first in Phase 0 (2026-11-21): the run never uses the
maintainer's checkout `~/src/SIMPLE` or its build and puts the machine's load beside every timing; no
speed verdict is taken on a loaded machine without the load beside it (R1). `/data/bgal`, the earlier
runs' directories and `~/tmp/blas_review` are read only; a bgal project is copied into the run's scratch
before a run writes it (R3). The base is the copy's HEAD, `7a679f424` as pinned above (R4).

Stops: decision 12's premise check; the kernel tolerance of decision 11 not met after diagnosis; a
bgal map outside the tolerance of section 6 (rulings R5) after an adoption; Phase 1 not byte-identical after diagnosis (the
shared tail must not change arithmetic); a file outside the table; the block order or the ordering of
decision 18 found to change a draw or a particle's candidates (a keyed draw depends on no order: if
one does, stop). A neighbourhood path failing its bar is a recorded negative result, not a stop.
Gates no ruling may waive: the fast gate after every edit, with the tester of decision 11 in it;
`suite=fills` identity; Phase 1's 6vxx identity; the 6vxx profiles of Phases 3 and 6.

## 8. Phases

- Phase 0, baseline (measurement, no edit).
  Done 2026-11-20 on the Dell at `f9fd7ca21` (Release, gfortran 15.2.1, bundled OpenBLAS 0.3.34
  OpenMP flavour, 8 P and 16 E cores): the benchmark table of section 1, the tile sweep, the check
  (`matrix_t` versus replica 4.0e-7 relative, best rotation 2 of 320,000, team identical to one
  thread in five repeats). The run's Phase 0 adds, at its base: the base build kept; the fast gate;
  `suite=fills` on the base; the bgal `solve3D` from `/data/bgal/tst/3_selection/tst.simple` twice
  with the same seed (the command of `phases.md`), keeping both runs' `BENCH/` records, maps and FSC;
  `test=profile target=simulated_workflow_6vxx nseeds=20 nthr=1 ncunits=24` as the before-profile;
  from the bgal records, per stage and table, the fill, its score loop and refinement, the references
  scored per particle, and the representatives' share of the `state` mode's pairs.
  Exit: the numbers above in Progress with their paths; the base's two-run spreads for the bgal
  gates; decision 12 checked.
- Phase 1, the kernel in the calculator.
  The shared tail (decision 14); the derived memo with its release (decision 6) and the memory term's
  routine (not yet called); the `cgemm` interface; the new submodule with `score_block_posterior`:
  the block's particle terms (memo or seed-bank rows), consecutive or listed tiles, the products, the
  transpose, the batched c2r, then per pair the shared conversion and draw with the pair's keyed
  variate; the per-thread work and the plan at `new`/`kill` (decision 7); the tester (decision 11);
  the harness's objective (decision 15).
  Exit: fast gate 15 of 15 with the new tester passing its checks inside one second;
  `unit_pftc_align2D3D` passes; `test=identity` of 6vxx at one thread against the base build
  identical; the harness's check re-run and recorded.
- Phase 2, the dense fill.
  `fill_tab_range` in blocks (decisions 4 and 5), the memo asked for in `new_pipeline`, the
  refinement reading `rows_batch`, the cost record's `ts_score` from the block score, the memory term
  counted when the dense table runs.
  Exit: `suite=fills` zero differing entries; the bgal run twice: decision 10 met (as restated by
  rulings R6), map and FSC within section 6's tolerance (rulings R5), the fill records recorded; fast
  gate.
- Phase 3, the dense profile.
  The 6vxx after-profile (20 seeds, one thread) against the before-profile.
  Exit: inside the spread, no gate row failing.
- Phase 4, the neighbourhood measured.
  The benchmark's neighbourhood route (decision 19), paths (a), (b), (c) on the kernel of Phase 1;
  their checks; decision 17 (i) per path.
  Exit: per path the pair counts, wasted pairs, thread-seconds against per-pair, the check inside
  decision 11's tolerances (the best rotation and the distances; the drawn rotations reported, rulings
  R7), and the verdict of decision 17 (i) in Progress; fast gate.
- Phase 5, the neighbourhood adopted.
  The paths that passed Phase 4 in `fill_tab_subspace_range` (decisions 16 and 18), the memo asked for
  in the neighbourhood table's `new_pipeline`, the memory term counted when it runs. If none passed:
  a Progress entry, no edit.
  Exit: `suite=fills` zero differing entries; the bgal run twice: decision 17 (ii) met per adopted
  path (as restated by rulings R6), map and FSC within section 6's tolerance (rulings R5); fast gate.
- Phase 6, the workflow profile.
  The 6vxx after-profile with both tables against Phase 0's before-profile; `multistate_workflows
  suite=labels` once. Skipped with a Progress entry when Phase 5 adopted nothing.
  Exit: inside the spread, no gate row failing; the labels route's accuracy reported.
- Phase 7, docs and report.
  `refine3D_policy.md` section 8 (the fills' blocks, the derived memo and its term, which
  neighbourhood paths use the kernel), a pointer from the benchmark report's "What it would take" to
  this note, the plan's status line and its README row, the report
  `prob_table_matrix_scoring_report.md` beside this note with the pftc skill's proposed lines and the
  commit message.
  Exit: the files of the table changed and nothing else; the report written.

### File table

| File | What | Phases |
|---|---|---|
| `src/main/pftc/simple_polarft_calc.f90` | components (memo, per-thread block work, plan) and interfaces, the shared tail's declaration | 1, 5 |
| `src/main/pftc/simple_polarft_block.f90` | new: `score_block_posterior` | 1, 5 |
| `src/main/pftc/simple_polarft_core.f90` | per-thread block work and plan at `new`/`kill` | 1, 5 |
| `src/main/pftc/simple_polarft_memo.f90` | the derived memo and its release | 1, 5 |
| `src/main/pftc/simple_polarft_corr.f90` | the shared tail as module procedures (decision 14); no arithmetic change | 1 |
| `src/main/pftc/simple_polarft_corr_tester.f90` | the kernel's tester (decision 11) | 1, 5 |
| `src/utils/math/simple_linalg.f90` | the `cgemm` interface | 1 |
| `src/main/commanders/test/simple_polar_matching_bench.f90` | the Nyquist bin (decision 15); the production kernel as a layout if needed; the neighbourhood route (decision 19); the fills harness's `state` and `geom` routes (added in Phase 5) | 1, 2, 4, 5 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | the benchmark's driver: the production kernel's layout and its check against the per-pair draws; the neighbourhood route's cases (added in Phase 1) | 1, 2, 4 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the benchmark's help text: the production kernel's layout (Phase 1), the neighbourhood route | 1, 4 |
| `src/main/strategies/search/simple_matcher_memory.f90` | the memo's term, counted when a kernel table runs | 1, 2, 5 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the memo asked for in `new_pipeline`; the dense fill in blocks; the cost record | 2 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | the dense table asks for the memo (`new_pipeline`), and the kernel's memory term is counted when it runs (added in Phase 2) | 2, 5 |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | the dense table's other production caller (`learn_pose_student`) asks for the kernel's memo (added in Phase 7: the report's review found it missed since Phase 2) | 7 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the adopted neighbourhood paths in blocks; the memo asked for | 5 |
| `doc/policies/3D/refine3D_policy.md` | section 8 | 7 |
| `doc/implementation_notes/completed/polar_matching_layout_bench_report.md` | pointer to this note | 7 |
| `doc/implementation_notes/README.md` | this note's row | 7 |
| `doc/implementation_notes/planned/prob_table_matrix_scoring.md` | this note: Progress | each |
| `doc/implementation_notes/planned/prob_table_matrix_scoring_report.md` | new: the run's report (beside this note) | 7 |

## Progress

- 2026-11-20, Dell, `f9fd7ca21`, before the run: the benchmark on the Dell (section 1) and its
  check, from `~/tmp/blas_review/pmb/polar_matching_bench.tsv` and `polar_matching_check.tsv`
  (`simple_test_exec test=polar_matching_bench suite=<per_pair|replica|matrix_t|check> nspace=2500
  nstates=2 nthr=<24|12|1>`; one thread on 200 particles). No code changed.
- 2026-11-20, rulings: Hans chose the recommended option for decisions 1, 3, 4, 6, 10, 11, 12 and 13.
  An architecture review before the rulings replaced a first draft's three parallel regions per
  batch by blocks as the unit of work (decision 4), a memo built inside `memoize_refs` with a second
  mirror rule by a derived memo (decision 6), and added decisions 13-15. It also found that the first
  draft's explanation of the 2.7e-3 difference between the harness and production (vectorization)
  was wrong: the harness drops the angular Nyquist bin that production sums since `9933db06c`.
- 2026-11-20, re-ruling of decision 1: the first draft covered the dense table alone on the
  statement that only it had a cost record. The bgal records show the neighbourhood table runs most
  of a single-state `solve3D` (stages 1-2 and 6-8) at a score cost like the dense table's. Hans ruled
  "dense, then neighbourhood", the neighbourhood adopted only past a bar fixed before its
  measurement: decisions 16-20, Phases 4-6. The note was renamed from `dense_table_matrix_scoring.md`
  and re-pinned to `7a679f424` (one documentation commit after `f9fd7ca21`, touching no file here).
- 2026-11-21, Phase 0, baseline (measurement; no source file changed). The Dell (Oracle Linux 8.10,
  Intel Core Ultra 9 285: 8 performance and 16 efficiency cores without hyperthreading, 62 GB), base
  `7a679f424`, the copy's HEAD and the plan's pin, so rulings R4 needed no re-pin. gfortran 15.2.1
  (gcc-toolset-15), CMake 3.26.5, the bundled OpenBLAS 0.3.34 (OpenMP flavour; `.deps/` copied read
  only from the maintainer's checkout into each source tree, so it was not rebuilt). Load 0.3 at the
  start; the load of every timing below is in its run's `load.log`, and nothing else heavy ran beside
  a timed job. Every tree is configured with `cmake -S <src> -B <tree> -DCMAKE_BUILD_TYPE=Release
  -DBUILD_TESTS=ON -DNICE=OFF -DUSE_ARCHOPT=ON -DFORTLS=OFF -DPython3_EXECUTABLE=/usr/bin/python3.12`
  (`-O3 -march=native`), then `make -j12`, `make install`, `scripts/run_fast_gate.sh <tree>`. Files
  changed: this note (section 7, the rulings at launch; this entry). Paths below are under
  `/home/elmlundho/agent_runs/table_gemm/scratch`; the scripts are in `keep/tools/`.
  - Base build: `keep/src_base` (`git archive 7a679f424` plus this note and its README row), built in
    `keep/build_base`, kept for every later comparison and never rebuilt. The run's tree `repo/build` is
    configured and built the same way. Fast gate 15 of 15 on both (`logs/build_base.log`,
    `logs/build_repo.log`).
  - Thread identity of the fills on the base (`simple_test_exec test=polar_matching_bench suite=fills
    nthr=24`, `keep/bench/p0_fills_base/`): the four probabilistic routes (the dense and neighbourhood
    3D tables, the dense and sparse 2D tables) at one and two states, each at one thread against the
    team of 24, a third of the batch, three parts, and three parts with the earlier parts' committed
    outputs perturbed before the later parts' fills: 48 comparisons, zero
    differing entries in each of the seven categories (distances, in-plane indices, shifts,
    temperatures, seeds, assignment, extras). Exit item met.
  - The bgal `solve3D` (beta-galactosidase, 4,967 particles, four parts of six threads, box 256, point
    group D2; the command of `phases.md`, whose 21 keys the base build accepts) twice with
    `SIMPLE_SEED=20261120` (`keep/bgal/p0_base_1`, `p0_base_2`: the per-iteration cost records
    `BENCH/`, final maps, FSC, manifest, statistics, log, load): wall 14:17 and 14:32. Both end at
    Fourier shell correlation (FSC) 0.143 at 3.75 Å and 0.5 at 4.29 Å. The stages: 1-2 the
    neighbourhood table in `shc` mode (iterations 1-40, 500 and 1,000 directions), 3-5 the dense table
    (41-75, 1,000 and 2,500), 6-8 the neighbourhood table in `state` mode (76-92, 2,500 and 5,000), as
    section 1 says. The table fills take 363 and 370 s of the run's wall (the slowest part per
    iteration): 284 and 291 s the dense table, 66 and 65 s `state` mode, 13 and 14 s `shc` mode.
  - **The two runs are byte-identical.** Final and post-processed maps, `fsc_state01.bin`, the sigma2
    store, every per-iteration resolution file, the particles' 3D orientations (`print_project_field
    oritype=ptcl3D`) and every count of the fill records agree bit for bit; only timings differ. A third
    run at five threads per part (`keep/bgal/p0_base_t5`, wall 15:15) is byte-identical as well: on this
    build `solve3D` reproduces itself across runs and thread counts at a fixed seed and part count. The
    two-run spreads of section 6 are therefore zero: 0 Å for the final FSC at 0.143, correlation 1 for
    the final maps. Exit item met as measured; but taken as written, the map gate of Phases 2 and 5
    becomes byte identity with the base, which a last-bits change cannot meet and section 2 rules out
    ("judged by profile, never by byte identity"). The run sets no other tolerance (section 6 is the
    maintainer's). For the ruling it records two perturbations of the base, run the same way:
    - three parts of eight threads, same seed (`keep/bgal/p0_base_p3`, wall 14:09): the parts' merge
      order changes the last bits; the run leaves the base's trajectory at iteration 4 and ends two
      iterations later (94). Final FSC 0.143 at 3.75 Å and 0.5 at 4.29 Å, as the base; correlation of
      its final map with the base's 0.9984 over the box (0.9982 inside the 180 Å mask; post-processed
      maps 0.9963); the two maps' FSC stays above 0.80 to Nyquist (2.55 Å).
    - another seed, 20261121, four parts of six threads (`keep/bgal/p0_base_s2`, wall 15:58): a
      different trajectory throughout. Final FSC 0.143 at 4.03 Å and 0.5 at 6.28 Å; its map lies in
      another symmetry frame of the same hand (correlation 0.35 with the base's as it is, 0.74 after
      `dock_volpair`, 0.47 docked from its mirror); docked, the two maps' FSC falls below 0.5 at 13 Å.
      A seed changes the solution far more than last bits do, so its spread would not gate a last-bits
      change.
    Without a ruling, Phase 2 stops at its map gate (section 7: a bgal map outside the base's spread).
  - Fill records per stage and table (`keep/bgal/p0_base_{1,2}_fill.tsv`, one row per iteration and
    part, parsed by record; medians over the stage's iterations of the four parts' mean, run 1 / run 2;
    a pair is one particle against one reference over every in-plane rotation):

    | stage | table | directions | shells to | iterations | fill wall (s) | score loop (thread-s) | refinement (thread-s) | references per particle | representatives' share | µs per pair and thread |
    |---|---|---|---|---|---|---|---|---|---|---|
    | 1 | neighbourhood, `shc` | 500 | 15 | 20 | 0.13 / 0.14 | 0.58 / 0.61 | 0 / 0 | 161 | - | 3.0 / 3.1 |
    | 2 | neighbourhood, `shc` | 1,000 | 16 | 20 | 0.46 / 0.51 | 2.32 / 2.51 | 0 / 0 | 407 | - | 4.6 / 4.9 |
    | 3 | dense | 1,000 | 26 | 11 | 4.21 / 4.61 | 7.83 / 8.38 | 17.0 / 19.0 | 1,000 | - | 6.3 / 6.7 |
    | 4 | dense | 1,000 | 33 | 15 | 4.47 / 4.65 | 9.88 / 10.28 | 16.0 / 17.3 | 1,000 | - | 8.0 / 8.3 |
    | 5 | dense | 2,500 | 30-65 | 9 | 15.1 / 15.3 | 75.3 / 76.3 | 18.7 / 18.6 | 2,500 | - | 24.2 / 24.6 |
    | 6 | neighbourhood, `state` | 2,500 | 65-70 | 6 | 2.15 / 2.13 | 4.36 / 4.29 | 7.01 / 6.94 | 186.5 | 68 % | 18.8 / 18.5 |
    | 7 | neighbourhood, `state` | 5,000 | 70-75 | 3 | 3.83 / 3.87 | 9.67 / 9.86 | 9.91 / 10.0 | 314.2 | 80 % | 24.8 / 25.3 |
    | 8 | neighbourhood, `state` | 5,000 | 75-81 | 8 | 4.50 / 4.36 | 11.7 / 11.1 | 11.8 / 11.7 | 314.2 | 80 % | 30.1 / 28.5 |

    The fills' load balance is 0.99 in the dense stages, 0.86-0.88 in `state` mode, 0.74-0.83 in
    `shc` mode. In `state` mode the subspace representatives are 126 per state at 2,500 directions and
    252 at 5,000 (the controller's `nint(126·nspace/2500)`), 68 % and 80 % of the references a
    particle scores; the other 60.5 and 62.2 are the members of about 3.1 cells of 19.8 references
    (the three peak cells, and the previous direction's cell when it is not among them), which agrees
    with those counts.
  - Decision 12 (stop before any edit if the dense score loop is below a quarter of the dense fill in
    every dense iteration): the score loop's share of the fill's thread-seconds per dense iteration
    (parts summed) is 0.18-0.94, median 0.36 (run 2: 0.17-0.94, median 0.34); one of the 35 dense
    iterations is below a quarter in both runs (iteration 41, the first, whose shift refinement costs
    twice that of the stage's later iterations), 4 of the 140 iteration-part records. No stop. The share is about a
    third in stages 3-4 (1,000 directions, shells to 26-33, where the shift refinement dominates) and
    0.8 in stage 5.
  - Decision 10's base, for Phase 2 (the bar itself is frozen): over the 35 dense iterations the median
    per-iteration score loop (mean of the four parts) is 9.54 and 10.15 thread-s; the median fill region
    wall of the slowest part 5.04 and 5.46 s; the per-iteration two-run spread of that wall has median
    0.35 s (4.7 %) and maximum 1.85 s (26 %, iteration 41, the first dense iteration), per iteration in
    `keep/bgal/p0_base_spread.tsv`. 26 of the 35 dense iterations are stages 3-4, scored on shells to
    26-33, about half the benchmark's 61 shells; the median is theirs.
  - 6vxx before-profile (`simple_test_exec test=profile target=simulated_workflow_6vxx nseeds=20 nthr=1
    ncunits=24 keep_outputs=no dir=keep/profiles/p0_base`, 28 min; `keep_outputs=no` deletes a passed
    execution's outputs as it ends, 5.5 GB each, and changes no metric): 14 of 20 executions pass,
    failure rate 0.30. Medians: map correlation with the truth 0.892 (spread `s` 0.029), FSC 0.143
    against the truth 4.30 Å (`s` 0.26), pose error of particle pairs 0.94°, fraction of pairs within
    5° 0.984. The six failures are the base's and are left (section 6): seed offset 3 at a truth
    correlation of 0.783 against the floor 0.80, its only failed row; offsets 4, 7, 11, 16 and 19 wrong
    maps (truth correlation 0.48-0.56, FSC 36-62 Å, pair pose error 60-85°). The nightly triage note
    records 2 of 20 missed maps on 2026-10-19.
  - What it means, in short: the base builds, gates and fills are as expected, the dense score loop is
    worth the kernel (a third of the dense fill at 1,000 directions, four fifths at 2,500; the dense
    fills are a third of the bgal run), and the representatives are four fifths of `state` mode's
    pairs at 5,000 directions. The bgal `solve3D` is now reproducible bit for bit, so the map tolerance
    of section 6 measured as planned is zero; Phase 2's map gate needs the maintainer's ruling (the
    perturbations above are the material for it) or Phase 2 stops there. Exit: the numbers in Progress
    with their paths, the base's two-run spreads (zero, with the two perturbations beside them) and
    decision 12 checked: met.
- 2026-11-21, Phase 1, the kernel in the calculator (identical; no production caller yet). The Dell, base
  `7a679f424`, the run's tree `repo/build` (configured as in Phase 0); load 1-5 outside the timed runs.
  Files changed: `simple_polarft_calc.f90`, `simple_polarft_corr.f90`, `simple_polarft_memo.f90`,
  `simple_polarft_core.f90`, new `simple_polarft_block.f90`, `simple_polarft_corr_tester.f90`
  (`src/main/pftc`), `simple_linalg.f90`, `simple_matcher_memory.f90`, the benchmark
  `simple_polar_matching_bench.f90` and, added to the file table before their edits, its driver
  `simple_commanders_test_highlevel_matchbench.f90` and its help text in `simple_test_ui_highlevel.f90`;
  this note (decisions 2 and 7, the file table, this entry). Rulings R2 are applied here for the first
  time (decisions 2 and 7). Paths below are under `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - The shared tail (decision 14): the conversion of a pair's inverse transform into distances,
    `euclid_dists_from_crvec`, is now a separate module procedure declared in the calculator's interface,
    and the draw is `draw_posterior_rvec` over an explicit real vector, `draw_posterior_crvec` its
    one-line caller on the thread's buffer. Rebuilt, fast gate 15 of 15 (`logs/rebuild_p1_tail.log`);
    the production kernels' distance vectors (the benchmark's `suite=identity`, 200 directions, 2 states,
    one thread) are byte-identical between the base build and this one, euclid and cc
    (`keep/bench/p1_kident_{base,run}`).
  - The kernel: `score_block_posterior` (new submodule `simple_polarft_block.f90`) gives, for up to 16
    particles of one half and any list of references, the entries `gen_posterior_val` (at zero shift) or
    `gen_posterior_val_bank` (per particle, from its seed bank) give: one `cgemm` per angular bin
    1..pftsz+1 per tile of 32 references, depth 2·nk, the products transposed per pair, one batched
    inverse FFT per particle's tile, then the shared tail with the pair's keyed variate. The bin-major
    memo `ft_ref_bins` is derived by `memoize_ref_bins` after `memoize_refs` (released with the
    reference memos and at `kill`; a generation counter stops the kernel on a stale memo); the per-thread
    work and the plan are made with the calculator's other plans and are not touched until used; an
    explicit `cgemm` interface is public in `simple_linalg`; the memory term `matcher_block_kernel_bytes`
    is defined and not yet called.
  - Why the product's shape is fixed (decisions 2 and 7 as revised): probes against the bundled OpenBLAS
    (`keep/blas_probe/probe*.f90`) show that `cgemm` gives an element other last bits when it sits in
    another column or the call has another shape (72 % of the elements of a 32 by 16 product differ from
    one-by-one products; 28 % when rows and columns are permuted in a full product), but not when only
    its row moves in a full product. With the particles on 16 fixed rows (absent ones zero) and the
    references on 32 fixed columns at their place in the list, 0 of about 160,000 pairs differ when
    their rows move, when padding changes or when the other columns' contents change, at one and 24
    threads. Without this a pair's last bits would depend on its block-mates, which decision 3 makes
    depend on the team and the batch.
  - The tester (decision 11), `test_block_kernel` in the polar correlation suite: 200 directions at one
    and two states, 32 particles of both halves with astigmatic CTFs, with and without seed banks, with
    consecutive and listed (shuffled) tiles, against `gen_posterior_val` and `gen_posterior_val_bank`:
    the drawn rotation differs in 0 of 6,400 or 12,800 pairs in each of the eight cases (allowed 1), the
    drawn loss within 2.5e-7 relative of the per-pair loss (bound 1e-6); with consecutive and with listed
    tiles a pair's entry is the same for blocks of 6, of 16 and of one, and for chunks of 64 references,
    and a forced team (four in the gate) gives the one-thread entries bit for bit in three repeats. Cost
    0.47-0.49 s at one thread (three runs, `keep/tester/t1b_rep*.log`; budget 1 s). Fast gate 15 of 15
    (`logs/rebuild_p1_review.log`). Stress by hand at a team of 24, ten runs of three repeats each, both
    tile kinds: identical every time (`keep/tester/t24b_rep*.log`).
  - The benchmark (decision 15): its replica and matrix layouts now sum the angular Nyquist bin. Check at
    the base sizes of 2026-11-20 (`suite=check nspace=2500 nstates=2 nthr=24`, `keep/bench/p1_check2`):
    the replica's largest difference from production fell from 2.66e-3 to 1.2e-7 (one unit in the last
    place) and the matrix layouts' from 2.67e-3 to 6.0e-7; `matrix_t` against the replica 4.0e-7
    relative with the best rotation differing in 2 of 320,000 pairs, as before; every layout
    thread-identical. A new layout `block_prod` runs production's block kernel as the dense fill will;
    against production's per-pair draws on the same keyed variates (64 particles, 5,000 references): the
    drawn loss within 3.3e-7 relative, the entries identical at 24 threads (five repeats) and for blocks
    of 4, 8 and 16; the drawn rotation differs in 70 of 320,000 pairs at the particles' own
    temperatures (411-429), in 5 at a tenth of them, in 2 at a huge temperature (the argmin). The draw
    multiplies the last-bit differences of the losses by the temperature: the argmin agrees at the rate
    decision 11 was measured on (2 of 320,000), the draw at these temperatures does not. Decision 11's
    tester is not affected (its draws: 0 differing), but decision 19's checks, held to decision 11's
    tolerance on the draw at the late-stage sizes, will meet the draw's rate; Phase 4 should compare
    the argmin, or the maintainer rules.
  - Speed at the base sizes (informational, `keep/bench/p1_timing/timing.log`, the machine otherwise
    idle; load as recorded per row), µs per pair and thread: `per_pair` 29.6 at 24 threads, 26.9 at 12,
    15.8 at one; `matrix_t` (B 16, R 32) 9.7, 7.7, 7.5; `block_prod` (B 16) 10.5, 7.6, 6.9, that is 2.8,
    3.5 and 2.3 times `per_pair`; at B 8 15.0 and at B 4 20.6 at 24 threads, since the fixed product
    shape computes 16 rows whatever B.
  - Identity: one execution of `simulated_workflow_6vxx` per build at one thread and seed offset 0
    (`test=profile ... nseeds=1 nthr=1`; the base build `keep/build_base`, 795 s; this phase's frozen
    build `keep/frozen/p1_kernel`, 713 s; both passed), then `test=identity`: the 185 declared outputs
    identical, 0 differing or missing (`keep/identity/cmp/identity_p1.log`).
  - Exit: fast gate 15 of 15 with the new tester passing its checks in under one second, met;
    `unit_pftc_align2D3D` passes (in the gate), met; the 6vxx identity at one thread against the base
    build, met; the harness's check re-run at the base sizes and recorded above, met. One comment of the
    new submodule changed after the identity run (no code); rebuilt, fast gate 15 of 15
    (`logs/rebuild_p1_final.log`).
- 2026-11-21, Phase 2, the dense fill in blocks: built, thread-identical, and stopped at two of its exit
  items (the map against the base's zero spread, and decision 10; both answered by rulings R5 and R6,
  see the next entry, which supersedes this entry's verdicts). The Dell, base `7a679f424`; the
  phase's frozen build `keep/frozen/p2_dense` (`diff.patch` beside it). Files changed:
  `simple_eul_prob_tab.f90`, `simple_strategy3D_matcher.f90` (added to the file table for Phases 2 and 5
  before its edit), the benchmark `simple_polar_matching_bench.f90` and its driver, this note (the file
  table, this entry). Paths below are under `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - The fill (decisions 3, 4, 5): `fill_tab_range` orders a batch's particles by half and cuts each half
    into blocks of min(16, max(4, n/(2·team))); one parallel region per batch, dynamic over blocks; a
    block is processed whole by one thread: each particle's seed (with the seeded path its seed bank, in
    the block's slot of `seed_bank(:,:,16,team)`), the block's rows scored at once by
    `score_block_posterior` on every table reference, then each particle's remaining steps in its path's
    order (the refinement, the shift map's screen and refinement, the temperature, the diagnostic). The
    per-thread copies of a row are gone: the refinement and the screen read the batch's rows. The search
    object's shift prior is set again before a particle's refinement (the block's later seeds set
    theirs). `new_pipeline` asks the calculator for the bin-major memo when called with `block_memo`
    (the matcher's dense table and the fills harness do; the neighbourhood table and the testers do not);
    the matcher counts the kernel's memory term when the dense table runs. The cost record's score
    thread-seconds are the block kernel's.
  - Thread identity: `suite=fills` (`keep/bench/p2_fills`, `p2_fills6`, `p2_fills6_rep2`): zero
    differing entries in every comparison, including two routes added to the harness because its dense
    route runs with the polar shift map only: `3D_dense_seeded` (`shift_sigma=0`, the seed and its
    bank) and `3D_dense_noshift` (`trs=0`); 72 comparisons per run (one thread against 24, a third of the
    batch, three parts, perturbed parts), the 24-thread configurations repeated three times. Fast gate
    15 of 15 (`logs/rebuild_p2_routes.log`).
  - The bgal `solve3D` twice on the frozen build, a base run between them (`keep/bgal/p2_*`, per-run
    TSVs and summaries beside them, the verdict in `keep/bgal/p2_decision10_verdict.txt`): walls 12:57
    and 13:18 against the base's 14:17, 14:32 and 14:33 (load from 2 to 13 at the starts, 15-18 during,
    the runs' own). The two runs of the new build are byte-identical to each other, as the base's are;
    the base run between them is byte-identical to Phase 0's. Final FSC at 0.143 3.75 Å and at 0.5
    4.29 Å, as the base. Per dense stage (run 1 / run 2 against the base's run 1): score loop
    thread-seconds stage 3 (1,000 directions, shells to 26) 4.42 / 4.88 against 7.83; stage 4 (shells to
    33) 4.93 / 5.39 against 9.88; stage 5 (2,500 directions, shells to 65) 23.4 / 23.6 against 75.3, that
    is 3.2 times faster where the sums are long. The table fill per iteration (mean of the parts): 3.84 /
    4.20, 3.74 / 3.94 and 7.33 / 7.56 s against 4.21, 4.47 and 15.1; the dense fills' sum over the run 179
    and 191 s against 259; their load balance 0.90 against 0.99 (a block is a coarser task than a
    particle).
  - Decision 10: the median per-iteration score loop over the 35 dense iterations is 4.92 and 5.37
    thread-seconds against the base's 9.54 and 10.15 (bar 4.92): 0.500 and 0.546 of the base, so run 1
    meets the bar to the last digit and run 2 does not; with the base run between them in the base the
    two are 0.48-0.50 and 0.52-0.55. No iteration's fill wall exceeds the base's by 5 % plus its two-run
    spread (the closest iteration 8.8 % and 16.8 % of the base below its limit). The median is set by stages 3-4 (26 of
    the 35 dense iterations), where the sums run over 25-32 shells, half the benchmark's 61; there the
    per-pair tail (the distances and the draw, unchanged and out of this plan's scope) and the inverse
    transform keep more of the cost, and decision 3's blocks of 12 particles (batches of about 150 per
    half on six threads) leave a quarter of the fixed 16-row product empty. Not met as written.
  - The map: the new build's final map correlates 0.99897 with the base's over the box (0.99886 inside
    the mask) and their FSC stays above 0.77 to Nyquist (`keep/bgal/p2_dense_1/cross_fsc_vs_base.log`).
    The base's own last-bits perturbation (three parts instead of four, Phase 0) correlates 0.99836 and
    keeps its FSC above 0.80. Against the base's two-run spread, which is zero because the base
    reproduces itself bit for bit, the map is outside: the stop of section 7.
  - What it means, in short: the dense fill in blocks is correct and thread-, batch- and part-identical,
    makes the dense fills 1.4 times faster over the bgal run and the whole run 9 % faster, and gives
    the base's resolution; its score loop is 3.2 times faster at the benchmark's sizes but only about
    two times faster at the short sums of stages 3-4, which set decision 10's median. Phase 2 stops for
    the maintainer's rulings on the map tolerance and on decision 10; the working tree keeps the change.
- 2026-11-21, Phase 2 resumed under rulings R5 and R6 (Hans, after the stop: "FSC unchanged within one
  shell, and correlation with the base at least 0.98 is fine"), applied here for the first time and
  recorded in decision 10, decision 17 (ii) and sections 5-8. No code changed since the stop: the tree
  is the frozen build's (`keep/frozen/p2_dense/diff.patch` and `status.txt` compared byte for byte), so
  the two bgal runs and the base runs of the previous entry are judged again, not rerun. Fast gate 15 of
  15 on the tree (`logs/gate_p2_resume.log`).
  - Decision 10 as restated (the dense score loop summed over every dense iteration and part, at most
    half of the base runs' mean sum): base 3,364 and 3,453 thread-seconds (mean 3,409; bar 1,704); the
    two runs 1,301 and 1,360, that is 0.38 and 0.40 of the base: met. With the interleaved base run in
    the base, 0.377-0.399. The wall-time guard as before: no iteration over in either run
    (`keep/bgal/p2_decision10_verdict_R6.txt`, by `keep/tools/dense_verdict.py`). For the record, the
    former per-iteration median stays at 0.500 and 0.546.
  - The map as ruled: both runs end at FSC 0.143 3.75 Å and 0.5 4.29 Å, the base's own values (zero
    shells apart); the final map, the same in both runs bit for bit, correlates 0.99897 with the base's
    over the whole map and 0.99886 within the mask (floor 0.98): met.
  - Exit: `suite=fills` zero differing entries, met (three runs: 48 comparisons on the four routes, then
    72 twice on six); the bgal run
    twice, decision 10 as restated met, map and FSC within the ruled tolerance, the fill records
    recorded (the previous entry and `keep/bgal/p2_*`), met; fast gate, met.
  - What it means, in short: the dense table now scores with the block kernel; over the bgal run its
    score loop costs 0.4 of what it did (0.31 at stage 5's long sums, 0.5-0.6 at stages 3-4's short
    ones), its fills take 1.4 times less wall time and the whole run 9 % less, with the base's
    resolution and a final map that correlates 0.999 with the base's.
- 2026-11-21, Phase 3, the dense profile. The Dell, base `7a679f424`; Phase 2's frozen build
  `keep/frozen/p2_dense` (the tree unchanged since, its tracked diff compared byte for byte). Files
  changed: this entry only. Paths below are under `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - The after-profile, run as the before-profile was (`simple_test_exec test=profile
    target=simulated_workflow_6vxx nseeds=20 nthr=1 ncunits=24 keep_outputs=no`, `keep/profiles/p3_dense`,
    1,704 s against 1,689; load 0.4 at the start, its own executions after). 14 of 20 executions pass,
    failure rate 0.30, the before-profile's: the same six seed offsets fail (3, 4, 7, 11, 16 and 19), each
    with metrics byte-identical to the before-profile's, the base's known failures (Phase 0).
  - The comparison (`test=profile dir1=keep/profiles/p0_base dir2=keep/profiles/p3_dense`,
    `keep/profiles/cmp_p0_p3/comparison_profile_simulated_workflow_6vxx_t1.tsv`): every gate row's median
    within its spread, shifted by 0 spreads except the map's correlation with the truth, +0.0002 spreads
    (0.89218 to 0.89219); of the report rows only the odd half's pair pose error moved, by -0.011
    spreads. 17 of the 20 executions write metrics identical to the before-profile's: at one thread the
    dense table runs 45 iterations of each `solve3D`, and where no drawn rotation or reference flips on
    the last bits the shift refinement after the draw returns the same poses. The three that differ
    (offsets 8, 14 and 15) improve slightly (truth correlation +0.0001, +0.0007 and +0.0051; offset 15's
    FSC against the truth from 4.38 to 4.30 Å). Median execution wall 763 s against 803.
  - Exit: inside the spread, no gate row failing beyond the before-profile's failure rate, met.
- 2026-11-21, Phase 4, the neighbourhood measured. The Dell, base `7a679f424`. Files changed: the
  benchmark (`simple_polar_matching_bench.f90`: the neighbourhood route, `polar_neigh_bench`), its driver
  (`simple_commanders_test_highlevel_matchbench.f90`: `suite=neigh`), its help
  (`simple_test_ui_highlevel.f90`), this note (decision 19, section 5, Phase 4's exit, this entry); no
  production source. Rulings R7 applied for the first time and recorded in decision 19, section 5 and
  the exit; the route's construction is recorded in decision 19. Fast gate 15 of 15 after each edit
  (`logs/rebuild_p4_neigh.log`, `rebuild_p4_ui.log`, `rebuild_p4_order.log`). Paths below are under
  `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - The runs: `simple_test_exec test=polar_matching_bench suite=neigh nthr=24` twice on frozen builds
    (`keep/frozen/p4_neigh`, then `p4_neigh_v2`, which adds a reported pass and scores otherwise the
    same), 228 and 350 s; per case (C1 and D2, one and two states) and path, two repeats at 24 threads
    and two at one thread, a team of four once, then the check. Load 1.5 at the start (a build's tail),
    then the run's own 2.3-7.4; nothing else ran. Evidence `keep/bench/p4_neigh` and
    `keep/bench/p4_neigh_run2` (the TSVs, the log, `verdict.txt` by `keep/tools/neigh_verdict.py`).
  - Pairs per particle (1,200 particles): (a) the representatives, 252 per state; (b) the members of the
    pooled cells, 60 at one state and 120 at two (bgal's stage 8 scores 315 references per particle,
    252 + 63); (c) the references inside the 10° cone, 39 and 77 at C1, 153 and 306 at D2. Wasted
    pairs, path (c) only: at C1 104,706 of 150,997 scored (2.3 times the needed) and 210,426 of
    303,022; at D2 194,787 of 378,195 (1.1 times) and 391,728 of 759,084. The kernel calls: (a) 100
    of 12 particles by 252 or 504 references, 74 % of the products' pairs used; (b) about 1,500 of 2.3-2.4
    particles by 20 or 40 references, 9-10 %; (c) 2,100-4,400 of 3.6-4.4 particles, 14-19 %.
  - Decision 17 (i), the block path's score thread-seconds over the per-pair routine's on the same
    pairs, summed over the two 24-thread repeats (the rule written in the journal before the runs),
    run 1 / run 2:

    | path | C1, 1 state | D2, 1 state | C1, 2 states | D2, 2 states | one thread | verdict |
    |---|---|---|---|---|---|---|
    | (a) representatives | 0.470 / 0.484 | 0.459 / 0.473 | 0.358 / 0.355 | 0.356 / 0.324 | 0.41-0.46 | passes |
    | (b) pooled members | 1.182 / 1.201 | 1.222 / 1.224 | 1.029 / 1.023 | 0.988 / 0.983 | 2.1-2.8 | fails |
    | (c) geom cones | 2.644 / 2.659 | 1.477 / 1.484 | 2.292 / 2.284 | 1.197 / 1.201 | 2.4-6.0 | fails |

    Per needed pair and thread at 24 threads (run 2): (a) 20-29 µs per pair, 9.5-9.9 in blocks; (b)
    52-53 µs per pair, 52-64 in blocks; (c) 52-53 µs per pair, 63-142 in blocks. Paths (b) and (c)
    are negative results (decision 20) and stay on the per-pair routine.
  - The check (rulings R7): best rotations differ in 1, 2, 9 and 2 of path (a)'s 302,400-604,800 pairs
    (allowed 6, 6, 12, 12) and in none of (b)'s and (c)'s; the distances agree within 4.8e-7 to 7.1e-7
    relative (tolerance 1e-6): inside decision 11's tolerances in every case, both runs. Reported: the
    drawn rotations at the particles' temperatures differ in 16-46 of (a)'s pairs and in none of (b)'s
    and (c)'s. Every pass, per pair and in blocks, is bit-identical to the first 24-thread pass at the
    second repeat, at one thread and at a team of four. The pools formed from (a)'s entries equal the
    per-pair routine's for all 1,200 particles in every case.
  - Why (b) and (c) fail: at 24 threads the per-pair routine scores a particle's members at 52 µs a pair
    against 14 at one thread, bound by memory (each particle reads its own 40-300 references); the
    blocks share little of it, because at the matcher's batch a cell is pooled by 2.4 particles of a
    block (reached by 3.6-4.4 cones), and the fixed 16-by-32 product (decision 2's revision) then uses
    a tenth to a fifth of its pairs, its particle terms rebuilt at every call. The cones add their
    outside pairs. bgal's parts batch 300 particles, so their blocks are sparser still.
  - Reported, not judged: the per-pair routine run over decision 18's blocks (one thread per block of
    particles ordered by previous cell) costs 0.86-0.93 of its run in particle order on the members and
    cones at 24 threads and 0.65-0.83 on the representatives; at one thread 0.98-1.08.
  - Projection for Phase 5, not a verdict (`keep/bench/p4_neigh_run2/projection_17ii.txt`): with (a) in
    blocks and the members per pair, `state` mode's score loop on this route is 0.68 of today's at one
    state (0.62-0.63 with the members scored in block order) and 0.53-0.56 at two; decision 17 (ii)'s
    bar of 2/3 on bgal's single-state run lies inside that range, so Phase 5's bgal runs decide it.
  - Exit: per path the pair counts, wasted pairs, thread-seconds against per-pair and the verdict of
    decision 17 (i), met (above: (a) passes, (b) and (c) fail); the checks inside decision 11's
    tolerances as rulings R7 restate them, met; fast gate, met.
- 2026-11-21, Phase 5, the neighbourhood adopted: path (a), the only one Phase 4 passed; (b) and (c)
  stay on the per-pair routine (decision 20). The Dell, base `7a679f424`. Files changed: the
  neighbourhood table (`simple_eul_prob_tab_neigh.f90`), the matcher (`simple_strategy3D_matcher.f90`),
  the benchmark's fills harness (`simple_polar_matching_bench.f90`, its file-table row extended to
  this phase before the edit), this note (decision 18's revision, this entry). Paths below are under
  `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - What changed: in `state` mode (`prob_neigh_mode=state`, stages 6-8 of a single-state `solve3D`)
    `fill_tab_subspace_range` now works in blocks. The batch's previous directions are found in one
    parallel pass; the particles are ordered by half and then by the subspace cell of their previous
    direction (decision 18 revised, recorded there) and cut into decision 3's blocks; one parallel
    region per batch runs dynamically over the blocks. A block is processed whole by one thread: each
    particle's shift seed and seed bank (in the block's slot), then every active state's subspace
    representatives scored for the whole block in one call of the block kernel, then each particle in
    turn: the representatives' entries recorded and its peaks kept in subspace order, its pooled
    neighbourhood, the members of the pooled cells scored per pair as before, its shift refinement
    (with its own shift prior set again). The per-pair coarse pass (`find_peak_subspaces`) lost its
    last caller and was deleted. `geom` mode keeps its per-particle region. `new_neigh_pipeline` takes
    the dense table's `block_memo` option: the matcher, whose references are memoized, passes it and
    the `state` mode then asks for the kernel's memo; the matcher counts the kernel's memory term when
    the `state` table runs.
  - Gates along the way: the fast gate 15 of 15 after every edit (`logs/rebuild_p5_*.log`). A first
    version asked for the memo inside `new_neigh_pipeline` unconditionally; the gate's probability-table
    tester builds a `state` table on a calculator without memoized references (and never fills it) and
    failed, so the request became the caller's, as for the dense table. A read-only review of the change
    (the representatives, temperatures and banks each particle is scored with, the order of its records,
    the search objects' state across a block's seeds, the slots and offsets, the cost counters, the
    parallel region) found no defect; two of its notes were acted on (`geom` allocates one bank slot
    per thread, not sixteen; a comment).
  - `suite=fills` with two new routes, the neighbourhood table's `state` and `geom` modes on a subspace
    map built as the builder builds it: before the adoption 96 comparisons, 0 differing
    (`keep/bench/p5_fills_pre`); on the frozen build (`keep/frozen/p5_state`) at 24 threads, again at 24
    and at a team of four, 3 x 96 comparisons, 0 differing, with the same entries per route as before
    the adoption (`keep/bench/p5_fills_post`): every configuration of one thread, the team, batches of a
    third and three parts gives the same rows, selections and extras.
  - The bgal `solve3D` twice on the frozen build (`keep/tools/bgal_run.sh`, the seed and the four parts
    of six threads of Phases 0 and 2), 12:57 and 13:11 (Phase 2's 12:57 and 13:18, the base's 14:17 and
    14:32); load 1.4 at the start, the run's own after (`keep/bgal/p5_state_{1,2}`, `run_both.log`).
  - Decision 17 (ii) as rulings R6 restate it: the `state` iterations' score loop summed over
    iterations and parts is 309.4 and 311.0 thread-seconds against the base runs' mean of 582.4 (bar
    388.3): 0.531 and 0.534 of the base, met (`keep/bgal/p5_decision17ii_verdict.txt`). Both runs have
    18 `state` iterations, as Phase 2's (the base 17); against Phase 2's runs, where the `state` mode
    was still per pair, 0.52. The wall-time guard: no iteration over, the closest 17.2 % and 18.3 % of
    the base below its limit. Per stage (median per iteration and part, run 1 against Phase 2's run 1):
    score loop 4.45 to 2.59 thread-seconds at stage 6 (2,500 directions), 9.98 to 5.12 at stage 7, 11.21
    to 5.52 at stage 8 (28.7 to 14.2 µs per pair); fill wall 2.17 to 1.88, 3.92 to 3.07 and 4.42 to
    3.42 s; the `state` fills' wall summed 61.3 to 49.3 s; load balance 0.87 to 0.78 (blocks of 12 are
    coarser tasks). The whole run's wall is unchanged within its noise. Phase 4's route projected 0.62-0.68;
    on bgal each of the four processes reads its own copy of the references, so the per-pair
    representatives cost more and the blocks save more.
  - The map: the final map is byte-identical to Phase 2's run's (MD5 sums in
    `keep/bgal/p5_state_1/final_map_md5.txt`), and every statistic the log prints for every iteration
    is the same: on bgal the representatives' last-bit differences change no peak, pool or draw. So the
    final FSC is the base's, 3.75 Å at 0.143 and 4.29 Å at 0.5 (zero shells apart), and the map
    correlates 0.99897 with the base's over the whole map and 0.99886 within the mask (floor 0.98): met.
  - Memory: the peak resident memory of a part at stage 8 grows from 2.93 to 4.64 GiB, the kernel's
    memo of all 5,000 references (decision 6) of which the `state` mode reads its 252 representatives;
    the matcher's estimate grows by the same 1.8 GB.
  - Exit: `suite=fills` zero differing entries, met; the bgal run twice, decision 17 (ii) met for path
    (a), map and FSC within section 6's tolerance, met; fast gate, met.
- 2026-11-21, Phase 6, the workflow profile. The Dell, base `7a679f424`; Phase 5's frozen build
  `keep/frozen/p5_state` (both tables on the block kernel: the dense table, and the neighbourhood table's
  `state` mode for its representatives). Files changed: this note only (rulings R8 in decisions 10 and
  17, this entry). Paths below are under `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - Rulings R8, applied here for the first time and recorded in decisions 10 and 17: the wall-time
    guard became the table's fill region wall summed over its iterations and parts, at most 1.05 times
    the base's. Phase 2's runs: dense fills 715.5 and 766.3 s against the base's mean of 1,064.0 s
    (0.672 and 0.720), met. Phase 5's runs: `state` fills 197.2 and 197.9 s against 237.0 s (0.832 and
    0.835), met; their dense fills 0.672 and 0.717. For scale, Phase 2's runs, whose `state` mode was
    still per pair, 1.035 and 1.040, and the third base run 0.991 (`keep/bgal/p6_R8_*.txt`, by
    `keep/tools/dense_verdict.py`).
  - The after-profile with both tables, run as Phases 0 and 3 ran theirs (`test=profile
    target=simulated_workflow_6vxx nseeds=20 nthr=1 ncunits=24`, `keep/profiles/p6_both`, 1,743 s; load
    0.8 at the start, its own executions after). Each execution's `solve3D` runs 56 `state` iterations
    (and 40 of the `shc` neighbourhood search), so the adopted path is exercised. 14 of 20 executions
    pass, failure rate 0.30, the before-profile's: the same six seed offsets fail (3, 4, 7, 11, 16 and
    19), the base's known failures (Phase 0).
  - The comparison with Phase 0's before-profile (`keep/profiles/cmp_p0_p6/comparison_profile_simulated_workflow_6vxx_t1.tsv`):
    every gate row's median within its spread, shifted by 0 spreads except the map's correlation with
    the truth (+0.0002 spreads, 0.89218 to 0.89219); of the report rows only the odd half's pair pose
    error moved (-0.011 spreads). This is Phase 3's comparison to the digit: every execution writes
    metrics byte-identical to Phase 3's (20 of 20), so the neighbourhood's adoption changes no
    workflow metric, as it changed no bgal output in Phase 5.
  - Reported, not gated: the median execution wall is 809 s (Phase 0's 803, Phase 3's 763); between
    profiles of the same code the per-execution walls move by up to 150 s either way (20 one-thread
    executions share 24 cores), and against Phase 3 the per-seed difference has a median of +16 s, so
    the wall shows no change beyond that noise. The peak resident memory per execution grows from
    2.26 to 3.58 GiB: the kernel's memo of every reference in the `state` stages (decision 6), of which
    the `state` mode reads only the representatives.
  - The multi-state label-recovery route once on the frozen build (`simple_test_exec
    test=multistate_workflows suite=labels nthr=8`, `keep/labels/labels`, 2:31; it runs 10 dense and 10
    `state` iterations on three states): every true state resolved, accuracy 0.999, adjusted Rand
    index 0.997, smallest recall 0.997, largest population error 0.001, median rotation error 1.80°
    (ceiling 2.0) and shift error 0.058 px (ceiling 0.10), the shuffled-label control resolving no
    state; all of the route's floors met. Against its profile of 2026-11-05 (dell2, recorded in
    `multistate_posterior_state_draw.md`, the directory not on this machine): the 13 of 20 runs that
    resolved every state had accuracy 0.998-1.000, and this run lies in that group. For context the
    same route on the base build (`keep/labels/labels_base`, 2:57) prints the same 41 result lines and
    the same per-iteration table, digit for digit.
  - Exit: the after-profile inside the before-profile's spread, no gate row failing beyond its failure
    rate, met; the labels route's accuracy reported, done.
- 2026-11-21, Phase 7, docs and report, and one fix. Files changed: `doc/policies/3D/refine3D_policy.md`
  (section 8: a new subsection on scoring in blocks, the `state` mode's line and the fills' schedule
  paragraph), `doc/implementation_notes/completed/polar_matching_layout_bench_report.md` (a pointer at
  the head of "What it would take in production"), `doc/implementation_notes/README.md` (this note's
  row), `src/main/strategies/parallelization/simple_strategy_pose_student.f90` (added to the file table
  with this phase before its edit), this note (the status line, the file table, Phase 2's fills counts
  corrected to 48, 72 and 72, this entry), and the new report `prob_table_matrix_scoring_report.md`
  beside this note (the findings, the changes per phase, the decisions taken under delegation, the
  deviations, the evidence against the base, what is left open, the pftc skill's proposed lines, the
  tree's `git status --porcelain` and a commit message). Paths below are under
  `/home/elmlundho/agent_runs/table_gemm/scratch`.
  - A read-only fact-check of the report against this Progress section and of the policy text against
    the code found the report's figures mostly right (the corrections are applied) and one regression:
    the pose student's training (`learn_pose_student`) builds a dense table after memoizing its
    references but without `block_memo`, then fills it, and has stopped at its first fill with "stale or
    missing bin-major memo" since Phase 2, which had searched only the matcher's call. Reproduced on
    Phase 5's frozen build (`simple_test_exec test=pose_student nthr=8`, `keep/pose_student/pose_pre`).
    The call now passes `block_memo=.true.`; rebuilt, fast gate 15 of 15
    (`logs/rebuild_p7_pose.log`), frozen as `keep/frozen/p7_final`. `test=pose_student` on it passes every
    floor, as on the base build (`keep/pose_student/pose_post`, `pose_base`, 42 s each); nine report rows
    without floors differ by a few thousandths (the seeded fraction 0.113 against 0.117; the merged map's
    FSC at 12 Å against the truth 0.74192 against 0.74189). No other caller fills a dense or `state` table
    without the memo (the matcher and the fills harness pass it; the class student and the 2D matcher use
    the 2D table).
  - Fast gate 15 of 15 on the final tree (`logs/gate_p7_final.log`).
  - Exit: the files of the table changed and nothing else, met (the pose student's strategy and the
    report added to the table before their edits); the report written, met.
