# Multi-state prob in stages 3-5, the table profile, the polar benchmark on the Dell: report

Run `ms_prob5check`, Phase 7 of `multistate_solve3D_regression.md` (beside this report), on the Dell,
2026-10-22. The code measured is HEAD ee5b68be1, which holds the maintainer's commits of 2026-10-21/22: no stage
skipping for multi-state runs, multi-state `prob` in stages 3-5 with `prob_neigh` only from stage 6, and the
vectorized per-pair polar kernels. No production code changed. The only code change is on the test side: the
benchmark's middle thread team. The diff is `/home/elmlundho/agent_runs/ms_prob5check/review/phase_7.diff` (written when the phase ends), the
journal `/home/elmlundho/agent_runs/ms_prob5check/state/journal_phase_7.md`.

## The three answers

1. **Does multi-state `solve3D` now separate the two states? Sometimes: one run of four separated, one more
   partly.** On Cyril's simulated set, from the 2D project made by the anchor build (P2), the four runs ended at
   accuracy 0.965 (separated), 0.806 (partly separated, still gaining when stage 5 ran out of iterations), 0.548
   and 0.510 (failed). The earlier counts were 0 of 4 on P2 for the build that ran `prob_neigh` in stage 5 (1 of 5
   with its run on the other 2D project) and 2 of 4 for the anchor build of 2026-10-15. Counts of four cannot rank
   these three builds; they would only show a difference as large as 4 of 4 against 0 of 4.

   What changed is where a split can still happen. With `prob` in stage 5, a split that starts late in stage 4
   carries on through stage 5: the separated run went 0.58, 0.71, 0.80, 0.88, 0.92 and on to 0.965, the partial
   run from 0.54 to 0.81. A run whose minority state drains in stage 4 does not recover, and neither does one whose
   labels freeze after one state locked on in stage 3. Every run found the structure (final pose errors about 1°).

2. **Is scoring most of the table time, and is the table most of the stage? The table is nearly all of the stage;
   scoring is about two thirds of the table.** In stages 3, 4 and 5 the probability table takes 86-93 % of each
   stage's wall time (77-94 % of single iterations). Inside it, scoring the particle-reference pairs takes 67 %, 57 % and 68 % of the thread
   time (four runs pooled; 51-88 % per run and stage), and the shift refinement of the best references takes
   nearly all the rest (31 %, 41 %, 31 %; 11-48 % per run and stage). The refinement's share moves with the data:
   the number of references refined per particle follows the spread of the previous poses, from 1 to 138 here.

   So scoring is about half of a prob iteration (49-60 %). Scoring 3.4-4.1 times faster, as the matrix-product
   layout does on this machine, would make a prob iteration about 1.5-1.8 times faster; the shift refinement would
   then be the larger part of the table.

3. **How do the kernels and layouts behave on the Dell? The new kernel is slower here, and the Dell gains more
   from reusing references than the Mac.**
   - **The vectorized euclid kernel of ee5b68be1 is slower than the old one on the Dell**, where it was 1.3
     times faster on the Mac. It runs at 0.92 times the old speed on one performance core (15.6 against 14.4 µs
     per pair), 0.81 times on one efficiency core (29.2 against 23.7 µs) and 0.75-0.96 times with 12 or 24
     threads (0.80-0.85 in a repeat in alternating order, which rules out an order effect). The cc kernel is 1.18 times faster on a performance core and 0.92 times on an
     efficiency core. On x86 the old loop was already vectorized with 256-bit vectors, so the Mac's cause does
     not apply here.
   - **The full machine is partly limited by memory, a single core by arithmetic.** On one core every layout
     costs the same per pair from 1,000 to 10,000 references. With all 24 threads production's cost per pair
     and thread doubles and grows with the number of references (31 to 43 µs), and in the real runs it rose from
     44 to 52-55 µs when stage 5 scored 5,000 references instead of stage 4's 2,000 at the same resolution.
   - **Layouts.** Scoring each reference against a block of four particles (`blocked`) is 1.7-2.1 times faster
     than production at 24 threads (1.6-1.75 times faster than the same arithmetic without the blocking), where
     on the Mac it gained nothing. The matrix-product layout (`matrix_t`) is 3.4-4.1 times faster than
     production at 24 threads and 3.8-4.9 times at 12, against 3.7 on the Mac. Its last-bit differences match
     the Mac's (at most 4.0e-7 relative, the best rotation changed in 0-3 pairs per case).

## Terms

- **State**: one of the two structures; each particle carries a state label. **Stage**: one rung of `solve3D`'s
  ladder of refine3D runs; stages 1-2 use `shc` (stochastic hill climbing), stages 3-5 now use `prob` (a dense
  table of particle-reference scores; the state is the best over all states).
- **Accuracy**: fraction of particles whose final state is their true state, under the better matching of found
  to true states. One state for all scores 0.567, a random split about 0.51. **Verdict** (fixed by the plan):
  separated at 0.90 or more, partial 0.75-0.90, failed below 0.75. **ARI**: adjusted Rand index (0 random, 1
  perfect).
- **Lock-on** (from the `ms_stage3diag` report): at the end of stage 3 one found state's map is right (docked onto
  a true map at correlation 0.93 or more, its particles' pose error a few degrees) and its half-maps agree to the
  edge of the 132-pixel crop (5.6 Å), which promotes stage 4 from 7.0 to 6.0 Å.
- **Pose error**: median angle between the relative orientation of two particles of the same true state and its
  true value (pairs within the same found state), so it does not depend on the frame or the hand.
- **Table profile**: the fill timers of the dense probability table (`prob_fill_profile`), written per part and
  iteration to `BENCH/PROB_TAB_BENCH_PART<p>.txt`: thread-seconds of the shift seed, the score loop (every pair:
  sums, inverse FFT, distances, in-plane draw) and the shift refinement of the best references, and the fill
  region's wall time. The iteration's wall time and its pre-step (the whole table step: three part jobs of 8
  threads scheduled by the master, then the table reads and the assignment) come from
  `BENCH/REFINE3D_STRATEGY_BENCH.txt`.
- **µs per pair per thread**: wall time × threads / pairs; lower is better. **P-core, E-core**: the Dell's 8
  performance and 16 efficiency cores.

## The four runs

Four runs of HEAD on P2, one at a time on the otherwise idle machine, each from its own copy of the project:
`simple_exec prg=solve3D pgrp=c1 mskdiam=270 nparts=3 nthr=8 nstates=2 print_states=yes`, without `SIMPLE_SEED`.
Kept for each run in `/home/elmlundho/agent_runs/ms_prob5check/scratch/keep/runs/H5_<n>/`: the labels after every
stage and iteration, the log, the final maps and their docking onto the true maps, the project and stage maps
after stages 2, 3 and 4 with pose errors and docking, the whole `BENCH/` directory, `metrics.txt`, `profile.txt`.

### Outcome

| run | final accuracy | ARI | populations (true 17,000 / 13,000) | verdict | final maps docked onto the truth |
|---|---|---|---|---|---|
| H5_1 | 0.806 | 0.371 | 7,716 / 22,284 | partial | map 2 onto true 1 at 0.985; map 1 onto true 2 at 0.936 but onto true 1 at 0.84 |
| H5_2 | 0.548 | 0.002 | 5,242 / 24,758 | failed | map 2 onto both at 0.97 (a mixture); map 1 onto neither (0.57, 0.61) |
| H5_3 | 0.510 | -0.006 | 24,025 / 5,975 | failed | both maps onto both truths at 0.96-0.98 (mixtures, in opposite hands) |
| H5_4 | 0.965 | 0.863 | 13,424 / 16,576 | separated | map 2 onto true 1 at 0.995, map 1 onto true 2 at 0.988 (both mirrored) |

Every run found the structure: the final pose error is 0.9-1.2° in all four (the median over pairs of particles
of the same true state). No run logged an empty state or a relabelling, and none omitted a stage or took an early
exit.

### Per stage

Accuracy after the stage, iterations, minutes. Stage 2 ends by its 20-iteration budget or by the scan exit (more
than 95 % of the search space scanned).

| stage | H5_1 | H5_2 | H5_3 | H5_4 |
|---|---|---|---|---|
| 1 | 0.50, 12 it, 1 min | 0.51, 12 it, 1 min | 0.52, 12 it, 1 min | 0.51, 12 it, 2 min |
| 2 | 0.50, 13 it (scan exit), 9 min | 0.51, 16 it (scan exit), 11 min | 0.51, 20 it, 14 min | 0.50, 14 it (scan exit), 10 min |
| 3 | 0.52, 17 it, 25 min | 0.51, 17 it, 27 min | 0.51, 17 it, 26 min | 0.51, 17 it, 25 min |
| 4 | 0.54, 15 it, 32 min | 0.54, 15 it, 32 min | 0.50, 12 it, 33 min | 0.58, 15 it, 45 min |
| 5 | 0.81, 12 it, 91 min | 0.55, 12 it, 69 min | 0.51, 12 it, 66 min | 0.96, 12 it, 83 min |
| stage 4 limit | 7.0 Å | 7.0 Å | 6.0 Å (promoted) | 6.0 Å (promoted) |
| run | 160 min | 141 min | 140 min | 165 min |

After stages 2, 3 and 4 (pose error: median for true state 1 / true state 2; docked: the best correlation of a
stage map onto either true map):

| run | after stage 2 | after stage 3 | after stage 4 |
|---|---|---|---|
| H5_1 | 31° / 28°, docked ≤ 0.43 | 23° / 16°, ≤ 0.58 | 2.3° / 2.6°, map 2 0.96 on both (mixture) |
| H5_2 | 31° / 29°, ≤ 0.23 | 24° / 18°, ≤ 0.56 | 2.1° / 2.3°, map 2 0.96 on both (mixture) |
| H5_3 | 24° / 18°, ≤ 0.55 | 8.2° / 2.7°, map 1 0.97 (mixture) | 1.3° / 1.4°, both 0.94-0.98 |
| H5_4 | 27° / 24°, ≤ 0.50 | 22° / 15°, ≤ 0.70 | 1.3° / 1.4°, both 0.96-0.98 |

### How each run went

- **H5_1, partial.** No lock-on in stage 3. The poses became right during stage 4 and the states started to split
  in its last four iterations (state 1 from 12,700 to 8,600 particles). In stage 5 state 1 shrank to 3,900 and
  then regrew while accuracy rose from 0.55 to 0.81, by 0.026-0.029 per iteration at the end: the run ran out of
  stage-5 iterations while still separating, with true state 2 split between the two found states.
- **H5_2, failed.** Like H5_1 up to stage 4, but the minority state drained instead of splitting: 11,100 to
  7,800 particles in stage 4, 5,200 at the end, holding a mixture with poor poses (10-12° against 1° in the
  majority). Accuracy stayed at 0.54-0.55 through stage 5.
- **H5_3, failed.** The only run with a full lock-on in stage 3 (one map right at 0.97, poses 8° and 2.7°, the
  half-maps at the crop edge, stage 4 promoted to 6.0 Å). After a swing of about 1,000 particles early in stage 4
  the labels froze at 23,900 / 6,100 for the rest of stage 4 and all of stage 5; both final maps are mixtures. This is the earlier fixed build's fourth run (HF_4) again,
  without even its small gain in stage 5.
- **H5_4, separated.** Stage 3 ended with state 2's half-maps at the crop edge and stage 4 promoted, but without
  a right map (docked ≤ 0.70, poses 15-22°). The poses became right in stage 4, the populations swung (state 1
  15,400 to 20,500 and back to 16,800), and accuracy rose from 0.50 to 0.58 over its last eight iterations: the
  anchor's shape. Stage 5 completed the split (0.71, 0.80, 0.88, 0.92, then 0.94-0.965).

### Against the earlier runs

Separation counts on P2 with the same command (the earlier fixed build also had one run on the 2D project made by
the current build, P2H):

| build | what it runs in stage 5 | separated | partial | failed |
|---|---|---|---|---|
| anchor 1bbbdcc24 (2026-10-15) | `prob` | 2 of 4 (0.971, 0.931) | 1 (0.768) | 1 (0.512) |
| fixed f8195504b (no stage skipping) | `prob_neigh` | 0 of 4 on P2 (+1 of 1 on P2H) | 0 | 4 |
| HEAD ee5b68be1 (this run) | `prob` (dense, 2 × 2,500 directions) | 1 of 4 (0.965) | 1 (0.806) | 2 |

The lock-on of the `ms_stage3diag` report (one state right after stage 3 and stage 4 promoted to 6.0 Å) against
the outcome, all builds:

| | separated | partial | failed |
|---|---|---|---|
| lock-on after stage 3 | A_1, AS_1, EF_1, H5_4 (marker only) | AS_2 | HF_4, H5_3 |
| no lock-on | | H5_1 | AS_3, HF_1, HF_2, HF_3, H5_2 |

The new runs add two things to that picture. A run without a stage-3 lock-on can still split, late, when stage 5
uses `prob` (H5_1; HF_2 started the same in stage 5 under `prob_neigh` but had less room). And a lock-on is not
enough: in H5_3, as in HF_4, the minority's labels froze in stage 4. Stage 2 ran its full 20 iterations in H5_3
and stopped by the scan exit after 13-16 in the other three; with one separated run among the three short ones,
these runs do not support the earlier lead that the scan exit costs lock-ons.

What the counts can and cannot tell: 1 separated and 1 partial in 4 sits between the earlier fixed build's 0 of 4
on P2 and the anchor's 2 of 4. Counts of four separate only large differences (4 of 4 against 0 of 4); the three
builds cannot be ranked on these runs.

### Time

| | stage 3 | stage 4 | stage 5 | run |
|---|---|---|---|---|
| HEAD (H5_1-4) | 25-27 min | 32-45 min | 66-91 min (`prob`) | 140-165 min |
| fixed build (HF_1-4, EF_1) | 18-27 min | 29-37 min | 26 min (`prob_neigh`) | 83-100 min |
| anchor (A_1, AS_1-3) | 29-45 min | 54-82 min | 100-138 min | 206-272 min |

Stage 5 is where the time went: the dense table over 5,000 references (2,500 directions per state) costs 5.5-7.6
minutes per iteration against about 2.2 with `prob_neigh`, and stage 5 is now half of a run. Stages 3 and 4 cost
about what they did. The slower euclid kernel (below) makes scoring roughly 5-30 % slower at 12-24 threads on this
machine, up to a tenth of a prob iteration; a stage-3 iteration takes 88-95 s against 73-94 s in the fixed build, consistent with that
but within the spread between runs.

## The table profile

From each run's `BENCH/` records (`scratch/keep/tools/profile.py`, per run in `profile.txt`), summed over the
iterations of each stage. "Table / iteration" is the master's pre-step (scheduling the three table parts, reading
their tables, assigning) against the iteration's wall time; "slowest part's fill" is the fill alone in the part
that finished last. The last four columns split the fill region's thread time (three parts of 8 threads).

| run | stage | iterations | stage wall, min | table (pre-step) / iteration | slowest part's fill / iteration | score | refinement | seed | rest |
|---|---|---|---|---|---|---|---|---|---|
| H5_1 | 3 | 17 | 25 | 89 % | 86 % | 68 % | 31 % | 0 % | 1.6 % |
| H5_1 | 4 | 15 | 32 | 88 % | 86 % | 60 % | 39 % | 0 % | 1.5 % |
| H5_1 | 5 | 12 | 91 | 93 % | 92 % | 52 % | 46 % | 0 % | 1.3 % |
| H5_2 | 3 | 17 | 27 | 90 % | 87 % | 62 % | 37 % | 0 % | 1.6 % |
| H5_2 | 4 | 15 | 32 | 87 % | 85 % | 64 % | 35 % | 0 % | 1.4 % |
| H5_2 | 5 | 12 | 69 | 91 % | 89 % | 76 % | 23 % | 0 % | 1.3 % |
| H5_3 | 3 | 17 | 26 | 89 % | 86 % | 66 % | 33 % | 0 % | 1.6 % |
| H5_3 | 4 | 12 | 33 | 87 % | 84 % | 57 % | 42 % | 0 % | 1.4 % |
| H5_3 | 5 | 12 | 66 | 86 % | 84 % | 88 % | 11 % | 0 % | 1.3 % |
| H5_4 | 3 | 17 | 25 | 86 % | 82 % | 75 % | 24 % | 0 % | 1.6 % |
| H5_4 | 4 | 15 | 45 | 88 % | 86 % | 51 % | 48 % | 0 % | 1.4 % |
| H5_4 | 5 | 12 | 83 | 89 % | 87 % | 63 % | 36 % | 0 % | 1.4 % |
| all runs | 3 | 68 | 103 | 88 % | 85 % | 67 % | 31 % | 0 % | 1.6 % |
| all runs | 4 | 57 | 141 | 88 % | 85 % | 57 % | 41 % | 0 % | 1.4 % |
| all runs | 5 | 48 | 310 | 90 % | 88 % | 68 % | 31 % | 0 % | 1.3 % |

What it says:

- **The table is the stage.** Summed over each prob stage of each run, the table step is 86-93 % of the wall time
  (77-94 % of single iterations), the fill alone 82-92 %. The rest (refine3D workers, assembly, the reference
  model) is 7-14 %. Inside a part the threads are balanced (0.98-0.99 of the region's thread time is used); the
  three parts finish within 4-8 % of each other (mean over slowest 0.92-0.96); the master's own share of the step
  (table reads, assignment) is about 1 %.
- **Scoring is the larger part of the table, not all of it.** Pooled over the four runs: 67 % of the fill's
  thread time in stage 3, 57 % in stage 4, 68 % in stage 5; the shift refinement of the best references takes
  31 %, 41 % and 31 %. Per run and stage the score share ranges from 51 % to 88 %. The multi-state table has no
  shift seed, and the remainder of the fill (batch setup inside the region, waits) is 1.3-1.6 %.
- **The refinement's share follows the data.** The table refines, per particle, the best references up to a
  number taken from the spread of the previous iteration's poses (`calc_num2sample` with `prob_athres`), the
  largest over the states: from 1 to 138 references per particle in these runs. It is large while a state is
  forming or splitting (H5_1's stage 5: 87-138, refinement 46 %) and small once the labels freeze (H5_3's stage
  5: 2-48, refinement 11 %).
- **Scoring costs more per pair in a real run than in the benchmark, and grows with the references.** The score
  loop costs 29-31 µs per pair and thread in stage 3 (2,000 references, 132-pixel crop), 37-44 in stage 4 (2,000,
  160 pixels) and 48-55 in stage 5 (5,000, 160 pixels). At the same resolution (H5_3 and H5_4 ran stage 4 at
  6.0 Å, as stage 5) going from 2,000 to 5,000 references raised it from 44 to 52-55 µs: memory again.
- **What a faster scoring would buy.** Scoring is about 49-60 % of a prob iteration (its share of the fill times the
  fill's share of the iteration, pooled over the runs). With the matrix-product layout's 3.4-4.1 times on this
  machine, a prob iteration would take about 0.55-0.66 of its time, 1.5-1.8 times faster, and stage 5 would drop
  from 66-91 to roughly 36-53 minutes (per run the score share varies, so the range is wider). The shift refinement would then be the larger part of the table.

## The polar matching benchmark on the Dell

### The machine

Intel Core Ultra 9 285: 8 performance cores (P, cpu0-7, each with its own 3 MB L2) and 16 efficiency cores (E,
cpu8-23, four to an L2), no hyperthreading, 36 MB L3, 62 GB DDR5. Oracle Linux 8.10 with kernel 4.18, which does
not know the hybrid layout: an unpinned one-thread job lands on either kind of core, so the one-thread rows of the
main series are a mixture. The pinned series (`taskset`) gives clean one-thread numbers per core kind.

Builds: gfortran 15.2.1, `-O3 -march=native` (AVX2), system OpenBLAS 0.3.15 (`OPENBLAS_CORETYPE=Haswell`,
`OPENBLAS_NUM_THREADS=1`), FFTW 3.3.5. The Mac: gfortran 16.2, `-O3 -mtune=generic`, OpenBLAS 0.3.27, FFTW 3.3.11.

### Old kernel against new kernel (production `per_pair`, `per_pair_cc`)

µs per pair per thread (lower is better), old (f3ccb5fbd) → new (HEAD), speedup = old / new.

| team | euclid, Dell | euclid, Mac | cc, Dell | cc, Mac |
|---|---|---|---|---|
| 1 thread, P-core (pinned) | 14.4 → 15.6 (0.92x) | 18.7-19.0 → 13.9-14.0 (1.35x) | 17.7 → 15.0 (1.18x) | 17.5-17.8 → 10.4 (1.70x) |
| 1 thread, E-core (pinned) | 23.7 → 29.2 (0.81x) | | 20.6 → 22.3 (0.92x) | |
| 12 threads (Mac 5) | 23-31 → 30-36 (0.75-0.86x) | 19.9-20.2 → 15.0-15.1 (1.33x) | 31-35 → 28-33 (1.05-1.14x) | 18.8-20.0 → 11.7-12.3 (1.61x) |
| 24 threads (Mac 11) | 29-37 → 31-43 (0.84-0.96x) | 30.9-31.1 → 23.9-24.1 (1.29x) | 27-39 → 31-36 (0.79-1.13x) | 28.8-30.4 → 19.8-20.7 (1.47x) |
| 8 P-cores (pinned) | 20.2-26.4 → 19.4-26.5 (0.88-1.04x) | | | |
| 16 E-cores (pinned) | 24-36 → 31-41 (0.67-0.88x) | | | |

The main series ran the new build before the old one in every suite. A repeat in the order new, old, new, old on
the smallest and the largest case gave the same picture (24 threads: 35.1 and 32.1 against 28.4 and 25.3 µs on
1,000 references, 40.0 and 42.7 against 34.5 and 36.1 on 10,000), and the pinned one-thread numbers repeated within 0.2
µs. It is not an order effect.

The vectorization reports of the two versions of `simple_polarft_corr.f90`, compiled with the build's flags, show
why the Mac's gain does not carry over: on x86 the old euclid loop was already vectorized with 32-byte (AVX2)
vectors behind an aliasing check. The pointer component that cost the Mac its de-interleaving loads did not stop
gfortran's x86 back end. What makes the new version slower here was not investigated further.

### The layouts on the Dell (HEAD)

µs per pair per thread and the speedup against production `per_pair`, best tile.

| layout | 24 threads | 12 threads | 1 P-core | 1 E-core | Mac 11 threads | Mac 1 thread |
|---|---|---|---|---|---|---|
| per_pair (production) | 31-43 | 30-36 | 15.6 | 29.2 | 23.9-24.1 (new) | 13.9-14.0 (new) |
| replica | 31-36 (1.01-1.27x) | 27-32 | 13.8 | 22.0-22.3 | 24.5-25.5 | 13.4-13.6 |
| blocked (B 4, R 8) | 18.5-21.1 (1.69-2.06x) | 19-21 (1.41-1.86x) | 10.9 | 20.2 | 23.8-25.1 | 13.3-13.5 |
| matrix_t | 8.7-10.7 (3.40-4.12x) | 7.1-8.0 (3.77-4.91x) | 4.5 | 6.4 | 8.3-8.4 (3.7x on old per_pair) | 5.2-5.7 |

Best `matrix_t` tiles on the Dell: B 16 R 32 (1,000-2,000 references), B 32 R 64 (4,000), B 16 R 128
(2,500-10,000), against B 32 R 64 everywhere on the Mac. Throughput at the full team: `per_pair` 0.56-0.77 million
pairs per second on the Dell (24 threads), 0.46 on the Mac (11 threads, new kernel); `matrix_t` 2.2-2.8 million on
the Dell, 1.3 on the Mac.

### Memory or arithmetic

- **One core is limited by arithmetic.** Pinned to one core, every layout costs the same per pair from 1,000 to
  10,000 references (within 0.3 µs, 1.5 %).
- **The full machine is limited by memory, in part.** Production's per-pair cost per thread at 24 threads is about
  twice the one-core cost and grows with the number of references, from 31 µs at 1,000 to 43 µs at 10,000. Eight
  P-cores alone already pay 19-26 µs against 15.6 µs on one, again growing with the reference count.
- **Reusing a reference pays here, unlike on the Mac.** `blocked` does the replica's arithmetic but scores each
  reference against four particles while it is in cache. Against the replica it is 1.6-1.75 times faster at 24
  threads, but only 1.27 times on one P-core and 1.09 times on one E-core. On the Mac it gained nothing over the
  replica at any team.
- **The matrix-product layout gains more here.** `matrix_t` is 3.4-4.1 times faster than production at 24 threads
  and 3.8-4.9 times at 12, and about 3.5 times on one P-core. It both reuses references and does the arithmetic in
  a faster form.

### Identity (suite=check, HEAD, 24 threads)

Every layout's five (or twenty) 24-thread repeats are bit-identical to its own one-thread run, in all six cases.
The replica, `blocked` and `states_inner` give identical distance vectors; their 20-lowest lists differ from
production's in 2.2-2.3 % of pairs (Mac 3.8-4.0 %), by the last bit. `matrix` and `matrix_t` differ from the
replica by at most 4.0e-7 relative, and the best rotation of a pair changed in 0-3 pairs per case (Mac 4.05e-7,
0-3).

## What changed, phase by phase

One phase, Phase 7 (`/home/elmlundho/agent_runs/ms_prob5check/review/phase_7.diff`):

- `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90`: the benchmark's middle thread team is
  half the full team instead of 5, the Mac's performance cores. On the Mac (11 threads) it is still 5; here it is
  12. Test side only; the fast gate passed after the change.
- The plan: two file-table rows for phase 7, the Phase 7 entry in Progress.
- This report.

## Decisions taken under delegation

All are in the journal; those that shape the evidence are also in the plan's Progress.

- **Release builds** (`compile_clean.sh`), not `compile_debug.sh`, because everything here is timed and compared
  with the earlier Release runs. Both builds passed the fast gate (15 of 15).
- **The middle team** of the benchmark became half the full team (above), the smallest change that gives the
  requested 1, 12 and 24 threads.
- **A pinned series** was added to the benchmark. The kernel 4.18 does not know the hybrid layout, so the
  unpinned one-thread rows jumped between about 15 and 29 µs depending on the core the job landed on. The pinned
  series ran every suite at one thread on a P-core (cpu3) and an E-core (cpu20), and `per_pair` on the 8 P-cores
  and the 16 E-cores. Its old-build `matrix_t` runs were skipped, because the benchmark's own layouts do not call
  the changed file and the main series had measured them equal on both builds.
- **A repeat in alternating order** (new, old, new, old) on the smallest and largest case, because the main series
  always ran the new build first and its result (the new kernel slower) was unexpected.
- **A vectorization report** of the two versions of the kernel file, compiled once with the build's flags into
  scratch, to see whether the Mac's cause applies on x86.
- **Snapshots after stage 4** as well as after stages 2 and 3, because in the first run the poses became good
  during stage 4, after the stage-3 marker.
- **Inputs read in place** from the earlier runs (P2, its stack, the true labels and maps, read only); the tools
  copied into this run's `scratch/keep/tools` and adapted to this build.

## Deviations from the plan

- The plan's phases 0-5 were not carried out; this run is its Phase 7 as the run's instructions define it.
- The benchmark's teams: 1, 12 and 24 as asked, plus the pinned series and the repeat (above).
- The runs took 2.3-2.75 hours each, not the 1.5 hours of the earlier fixed-build runs, because stage 5 now builds
  the dense table over 5,000 references. Nothing was dropped for it.
- No run was repeated, no threshold moved, no floor or assertion touched.

## Left open for the maintainer

- **The kernel commit on x86.** On the Dell, the kind of machine production runs on, ee5b68be1 makes a euclid pair
  take 8-23 % longer on one core and 4-33 % longer with 12 or 24 threads (18-25 % in the alternating repeat). Undoing the euclid part on x86,
  or finding a form that is fast on both machines, is a production decision this run could not take. The cc
  change is a gain on P-cores (1.18 times) and a small loss on E-cores (0.92).
- **The matrix-product layout is worth building, for the Dell more than for the Mac.** Scoring is about half of a
  prob iteration (49-60 %) and the layout scores 3.4-4.9 times faster here. Estimated gain: prob iterations
  1.5-1.8 times faster. The shift refinement of the best references (11-48 % of the table) is the next target.
- **Stage 5 is now the most expensive stage by far**: 66-91 minutes against 26 with `prob_neigh`, so a run takes
  2.3-2.75 hours against the 30 minutes ruling 1 expects. With `prob` in stage 5 two of the four runs separated
  fully or partly, against none of four on P2 with `prob_neigh`; whether that difference is real needs more runs.
- **Runs.** Four runs per build cannot separate 1 in 4 from 2 in 4. On this set every build tried so far (anchor
  2 of 4, the earlier fixed build 1 of 5, HEAD 1 of 4 with 1 partial) separates at most half of its runs, and the
  failures share two signatures: the minority state drains, or its labels freeze after a stage-3 lock-on.
