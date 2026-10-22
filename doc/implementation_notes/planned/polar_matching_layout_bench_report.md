# Polar matching layouts for many references: the report

Report of the measure-only benchmark planned in `polar_matching_layout_bench.md` (same directory). Unattended run on
Hans's Mac (Apple M3 Pro, 5 performance and 6 efficiency cores, 18 GB), 2026-10-21. No production code was changed.

## The answer

- **Scoring can be about 3.5 times faster.** The fastest layout computes the two sums of every angular frequency
  bin as one complex matrix product (BLAS `cgemm`) over a tile of references and a block of particles. It then
  transposes the products so that each pair's inverse FFT reads contiguous memory. This layout is called
  `matrix_t` below. Against today's per-pair routine it is 3.3-3.8 times faster, at 1, 5 and 11 threads alike.
- **It is the same for one, two and four states.** On this machine the speedup does not depend on the number of
  states or directions: 3.7 times at 11 threads, 3.5 at 5, 3.3-3.6 at 1. The cost per pair is flat from 1,000 to
  10,000 references.
- **The learning: on the Mac, scoring is limited by arithmetic, not memory.** The plan expected memory bandwidth to
  be the limit. It is not. Today's routine costs the same per pair whether the references fit in cache or not.
  Reusing a reference across several particles (layout 2) therefore gains almost nothing. What helps is doing the
  same arithmetic in a faster form: the matrix product.
- **The cost in memory is one copy of the reference spectra in another order.** The bin-major copy is as large as
  the memoized reference spectra it replaces for scoring: 3.6 GB for 2,500 directions and 4 states (both halves),
  0.9 GB for 2,500 directions and one state. It takes under a second to build. If production keeps the present
  arrays as well (the shift refinement reads them), the reference memory doubles.
- **The cost in identity is the last bits.** `matrix_t` sums in a different order and in single precision.
  Its distances differ from today's by at most 4.1e-7 relative. The best in-plane rotation of a pair changed in
  at most 3 of 320,000 pairs. The lists of the 20 best rotations, which feed the in-plane draw, differ somewhere in
  almost every pair, by at most 6e-7.
- **Layouts 2 and 3 are byte-identical but slow.** Blocking particles (layout 2) and batching the inverse FFT
  over states (layout 3) keep every distance bit for bit. They are 1.2-1.4 times faster than production. That gain
  is not theirs: production's own sums, compiled in the benchmark module without any change of layout (the
  "replica"), are as fast.

## Terms

- **Pair**: one particle scored against one reference over all in-plane rotations. One pair is two sums over the
  61 polar shells and 182 angular bins, one inverse real FFT of length 364, the conversion to 364 distances and the
  in-plane draw among the 20 lowest.
- **µs/pair/thread**: wall time × threads / pairs. Lower is better. At 11 threads it includes the slower
  efficiency cores.
- **Speedup**: pairs per second of a layout divided by that of layout 1, same case and threads.
- **Case**: number of projection directions (1,000 or 2,500) times number of states (1, 2, 4).
- **B, R**: particles per block and references per tile.
- **Replica**: production's per-pair sums (`gen_euclid_crvec`, unshifted branch) and FFT, the same source,
  compiled in the benchmark module on a memo the benchmark rebuilds. It is the reference for byte identity.

## The layouts

| name | what it does |
|---|---|
| `per_pair` (layout 1) | production `gen_likelihood_val`, called as `fill_tab_range` calls it: one particle per thread against every reference |
| `replica` | layout 1's sums and FFT on the benchmark's memo, particle by particle |
| `blocked` (layout 2) | a thread takes B particles of one half; per tile of R references, each reference against every particle of the block; per-pair arithmetic unchanged |
| `states_inner` (layout 3) | per particle and direction, every state's sums, then one batched inverse FFT over the states |
| `matrix` (layout 4) | per block and tile, one `cgemm` per angular bin, depth 2 × 61 (both sums fused), then R × B inverse FFTs whose input is strided |
| `matrix_t` (layout 4, variant) | as `matrix`, but the products are transposed per tile so that each FFT reads contiguous input |

`matrix_t` was added during the run when the profile showed that the strided FFT input of `matrix` grew expensive
for larger tiles.

## Sizes

Taken from SIMPLE's own parameter derivation: box 180 at 2.055 Å, crop 160, mask diameter 270 Å (crop radius 58),
low-pass limit 6 Å. That gives `pftsz` 182, `nrots` 364, shells 2-62 (61 shells), and 20 in-plane samples
(`prob_athres` 10°). Synthetic spectra: uniform random complex values with an exponential radial fall-off, a
physical astigmatic CTF per particle (through `create_polar_absctfmats`), and a noise power that falls to a floor,
slightly different per half. Particles alternate between the halves. Objective: euclid.

## The table

µs/pair/thread, and the speedup against layout 1. Best tile shown for `blocked`, `matrix`, `matrix_t`. 2,000
particles per measurement at 5 and 11 threads, 200 at 1 thread. Each value is the best of three after a warm-up.

| directions | states | refs | threads | per_pair | replica | blocked | states_inner | matrix | matrix_t |
|---|---|---|---|---|---|---|---|---|---|
| 1000 | 1 | 1000 | 11 | 30.9 | 24.5, 1.26x | 25.1, 1.23x (B4 R8) | 23.8, 1.30x | 11.8, 2.62x (B8 R32) | 8.3, 3.74x (B32 R64) |
| 1000 | 1 | 1000 | 5 | 20.0 | 15.0, 1.33x | 14.8, 1.35x (B4 R8) | 14.7, 1.36x | 7.6, 2.61x (B8 R32) | 5.8, 3.45x (B32 R64) |
| 1000 | 1 | 1000 | 1 | 18.7 | 13.4, 1.39x | 13.5, 1.39x (B4 R8) | 13.4, 1.40x | 6.9, 2.69x (B8 R32) | 5.3, 3.54x (B32 R64) |
| 1000 | 2 | 2000 | 11 | 31.1 | 24.8, 1.25x | 23.8, 1.31x (B4 R8) | 23.6, 1.32x | 11.3, 2.75x (B8 R32) | 8.3, 3.76x (B32 R64) |
| 1000 | 2 | 2000 | 5 | 20.0 | 14.7, 1.36x | 14.3, 1.40x (B4 R8) | 14.6, 1.37x | 7.6, 2.65x (B8 R32) | 5.8, 3.48x (B32 R64) |
| 1000 | 2 | 2000 | 1 | 18.9 | 13.5, 1.41x | 13.3, 1.43x (B4 R8) | 13.3, 1.42x | 6.9, 2.74x (B8 R32) | 5.3, 3.60x (B32 R64) |
| 1000 | 4 | 4000 | 11 | 31.0 | 25.0, 1.24x | 24.0, 1.29x (B4 R8) | 23.4, 1.33x | 11.2, 2.76x (B8 R32) | 8.3, 3.75x (B32 R64) |
| 1000 | 4 | 4000 | 5 | 20.0 | 14.8, 1.35x | 14.2, 1.40x (B4 R8) | 14.4, 1.38x | 7.7, 2.60x (B8 R32) | 5.7, 3.50x (B32 R64) |
| 1000 | 4 | 4000 | 1 | 18.7 | 13.6, 1.37x | 13.3, 1.41x (B4 R8) | 13.3, 1.41x | 7.1, 2.66x (B8 R32) | 5.2, 3.57x (B32 R64) |
| 2500 | 1 | 2500 | 11 | 31.1 | 25.5, 1.22x | 23.8, 1.30x (B8 R8) | 22.8, 1.36x | 12.7, 2.45x (B16 R32) | 8.3, 3.73x (B32 R64) |
| 2500 | 1 | 2500 | 5 | 20.0 | 14.6, 1.37x | 14.2, 1.40x (B8 R8) | 14.3, 1.40x | 8.2, 2.42x (B16 R32) | 5.8, 3.47x (B32 R64) |
| 2500 | 1 | 2500 | 1 | 18.7 | 13.6, 1.37x | 13.3, 1.41x (B8 R8) | 13.3, 1.41x | 7.6, 2.45x (B16 R32) | 5.3, 3.56x (B32 R64) |
| 2500 | 2 | 5000 | 11 | 31.1 | 24.6, 1.26x | 24.1, 1.29x (B8 R8) | 22.9, 1.36x | 12.9, 2.40x (B16 R32) | 8.3, 3.74x (B32 R64) |
| 2500 | 2 | 5000 | 5 | 20.0 | 14.6, 1.37x | 14.4, 1.39x (B8 R8) | 14.8, 1.35x | 8.7, 2.30x (B16 R32) | 5.8, 3.45x (B32 R64) |
| 2500 | 2 | 5000 | 1 | 18.7 | 13.6, 1.38x | 13.3, 1.41x (B8 R8) | 13.8, 1.36x | 8.1, 2.32x (B16 R32) | 5.7, 3.26x (B32 R64) |
| 2500 | 4 | 10000 | 11 | 31.1 | 25.2, 1.24x | 23.9, 1.30x (B4 R8) | 23.2, 1.34x | 13.1, 2.38x (B16 R32) | 8.4, 3.72x (B32 R64) |
| 2500 | 4 | 10000 | 5 | 20.4 | 14.4, 1.41x | 14.3, 1.42x (B4 R8) | 14.6, 1.40x | 8.2, 2.48x (B16 R32) | 5.8, 3.50x (B32 R64) |
| 2500 | 4 | 10000 | 1 | 19.1 | 13.4, 1.43x | 13.3, 1.43x (B4 R8) | 13.7, 1.39x | 8.2, 2.32x (B16 R32) | 5.3, 3.63x (B32 R64) |

Throughput at 11 threads, and the effective bandwidth (pairs per second × 178.6 KB of reference spectra one pair
reads):

| layout | pairs per second | effective GB/s | largest case wall time (10,000 refs × 2,000 particles) |
|---|---|---|---|
| per_pair | 3.53-3.56e5 | 63 | 56.6 s |
| replica | 4.32-4.49e5 | 77-80 | 45.8 s |
| matrix | 8.40-9.79e5 | 150-175 | 23.8 s |
| matrix_t | 1.31-1.33e6 | 234-238 | 15.2 s |

The effective bandwidth of `matrix_t` exceeds what the memory can deliver, because a reference is read once per
block of 32 particles. For `per_pair` it is well below the memory's capacity: it is not waiting for memory.

Tiles. For `blocked`, larger blocks are slower (B 4: 1.30x, B 32: 1.20x at 11 threads). For `matrix`, B 8-16 and
R 32 were best; larger tiles lost to the strided FFT. For `matrix_t`, B 32 and R 64 were best everywhere at 11
threads. B 32 was the largest block swept, so larger blocks may be faster still.

## Memory

| layout | reference spectra | extra |
|---|---|---|
| per_pair, replica, blocked, states_inner | `ft_ref` + `ft_ref2`: 357 KB per reference (both halves); 3.6 GB at 10,000 references | none |
| matrix, matrix_t | bin-major copy (bins 1-182): 355 KB per reference; 3.6 GB at 10,000 references | per thread about 15 MB at B 32, R 64 (particle block, products, transposed products, FFT output) |

The memo totals the program records, particles included, are 0.71 GB (1,000 references) to 3.9 GB (10,000).
Building the bin-major copy (the scatter, one thread) took 0.06 s for 1,000 references and 0.93 s for 10,000,
next to the FFTs the memo needs anyway.

## Identity

The checks score 64 particles (32 in the 10,000-reference case) against every reference. Each layout runs once at
one thread and five times at 11 threads.

| comparison | result |
|---|---|
| production at 11 threads (×5) against production at 1 thread | identical, all cases |
| every layout at 11 threads (×5) against its own 1-thread run | identical bytes, all cases (the stress test) |
| `blocked` (B 4, 8, 16, 32) and `states_inner` against the replica | identical distance vectors, all pairs, all cases |
| replica against production (the 20 lowest distances and their rotations, which `gen_likelihood_val` returns) | Release build: 3.8-4.0 % of pairs differ, by at most 1.2e-7 (one unit in the last place); Debug build: identical |
| `matrix`, `matrix_t` against the replica | largest relative difference 4.05e-7, absolute 6e-7; the best rotation differs in 0-3 pairs per case (of up to 320,000); `matrix` and `matrix_t` give the same differences |

Why the replica is not bit-identical to production in Release: the source of the sums is the same, the memo is
rebuilt with production's FFT plans, and in the Debug build the two agree exactly. A replica loop compiled
without vectorization differed in 33 % of pairs instead of 4 %; a constant CTF in both paths still differed.
The last bit of a distance therefore depends on how the compiler vectorizes the loop in its compilation unit.
The consequence for production: a new layout keeps byte identity only if it calls the very same compiled kernel.
"Same sums in the same order" in another routine is not enough.

The comparison is made on the distance vectors before the in-plane draw. The draw takes random numbers from each
thread's own stream (gfortran's `random_number`), so a different loop order changes which draws a pair meets, in
every layout.

## Where the time goes

Thread-µs per pair, from timers in the replica and the matrix layouts (2,500 directions; replica at 2 states and
400 particles, `matrix_t` at 4 states and 2,000):

| layout, threads | sums (or particle block + `cgemm`) | inverse FFT (with transpose) | distances and draw | total |
|---|---|---|---|---|
| replica, 1 | 10.8 | 0.64 | 2.3 | 13.8 |
| replica, 11 | 19.2 | 1.18 | 2.9 | 23.6 |
| matrix_t B32 R64, 1 | 2.3 | 0.76 | 2.3 | 5.3 |
| matrix_t B32 R64, 11 | 4.0 | 1.29 | 2.8 | 8.4 |

The sums shrink 4.7 times with `cgemm`. The distances and the draw (364 conversions in double precision, the
selection of the 20 lowest, the exponentials of the draw) are the same code in every layout. They are now
44 % of `matrix_t` at one thread. Any further gain has to come from there: faster sums alone cannot push past
about 2.3 µs per pair at one thread.

## What it would take in production

- **Where the code goes.** The memoized arrays are private to `polarft_calc`, so the kernel belongs there:
  - a bin-major reference copy built by `memoize_refs` (`simple_polarft_memo.f90`, components and interfaces in
    `simple_polarft_calc.f90`);
  - a block kernel in `simple_polarft_corr.f90` that scores B particles against R references. It would take the
    weighted particle terms, `cgemm` per angular bin, transpose, run a batched c2r, and return distance vectors;
  - an explicit `cgemm` interface, as `simple_linalg.f90` has one for `sgemm`.
- **The table fill.** `fill_tab_range` (`simple_eul_prob_tab.f90`) would score blocks of same-half particles
  instead of one particle per thread. Each pair's draw stays the existing `sample_likelihood_dists` call.
  - The shift seed fits. Since 2026-10-15 a particle's shift is applied on the particle side, so a particle's
    block row is built from its shifted memo (`memoize_shifted_ptcl`) and the references stay fixed. The
    single-state seeded path works the same way.
  - The shift refinement of the best references (`refine_best_refs`) is untouched; it still reads `ft_ref`. Both
    arrays then stay in memory: the reference memo doubles, +3.6 GB at 10,000 references, +0.9 GB at 2,500.
- **BLAS threads.** `cgemm` runs inside the OpenMP team, so the BLAS must be single-threaded. The benchmark runs
  with `OPENBLAS_NUM_THREADS=1`; production would need the same guarantee (environment or `openblas_set_num_threads`).
  With it, the anaconda OpenBLAS 0.3.27 gave identical bytes at 11 threads in every repeat.
- **The in-plane draw order.** The table entries would come out in block order. Each pair still draws from its
  thread's random stream, as today, so the result is statistically the same but not bit-reproducible against an
  earlier run. The 20-lowest lists that feed the draw differ in the last bits (above). Whether that matters is a
  decision about reproducibility, not about accuracy: 4e-7 relative is far below the spread of the likelihood.
- **Precision.** Production accumulates the sums with double-precision weights (`sigma2_wk8`). `cgemm` works in
  single precision. `zgemm` (double) would keep more bits at roughly twice the cost of the matrix part; it was not
  measured.

## The measurement on real data that should come first

This benchmark times the scoring of the dense table alone. It cannot say what share of a real run's table time
the scoring is. The shift refinement of the best references (L-BFGS from a coarse grid) depends on the data and
may dominate.

The measurement that settles it is the table's timing profile in the BENCH output of one real two-state
`solve3D` run on the Dell. `fill_tab_range` already records it (`prob_fill_profile`: `ts_seed`, `ts_score`,
`ts_refine`, `wall_region`). If `ts_score` is most of the table time, a 3.5-times faster scoring is worth the
production work. If `ts_refine` dominates, it is not.

The Dell should also run this benchmark (`simple_test_exec test=polar_matching_bench`). On the Dell the plan
estimated about 70 µs per pair and thread, close to streaming the reference spectra from memory. If the Dell is
limited by memory and the Mac by arithmetic, layout 2 may help there and the matrix layout's gain may differ.
That is a measurement, not a guess this report can make.

## Observations about production

No defect was found. Two observations:

- Production's per-pair sums run 1.2-1.4 times slower than the same source compiled in the benchmark module.
  Production reads its arrays through the derived-type components of `polarft_calc` (`self%crvec1(ithr)%c`,
  `self%ft_ref(...)`), which may limit vectorization. This was not investigated further. It may be worth a look
  on its own, independent of any change of layout.
- The memoized arrays have no getter, so a test cannot read them. The benchmark rebuilds them with production's
  FFT geometry, which agrees exactly in the Debug build.

## Caveats

- One machine. All timings are from the Apple M3 Pro Release build (gfortran 16.2, `-O3 -mtune=generic`,
  OpenMP, OpenBLAS 0.3.27 from anaconda, FFTW 3.3.11). FFTW plans were made with `FFTW_ESTIMATE` (fixed seed);
  production uses `FFTW_MEASURE` with wisdom.
- One thread used 200 particles per measurement instead of 2,000, to fit the time budget. The per-pair cost
  does not depend on the reference count, so this changes nothing measurable.
- Background load. From about 18:00 a video call ran on the Mac (about one core). `per_pair` re-timed at 19:51
  was 4.9 % slower at 11 threads and 1-2 % slower at 5 and 1 thread than in the first series. The other
  layouts were timed after 17:50, so their speedups at 11 threads may be slightly understated. macOS load
  averages read 40-130 during the large cases with no paging (12 swap-ins in 20 s); they do not reflect CPU use
  here.
- Synthetic data. The scoring cost does not depend on the values. The identity statistics (how often the best
  rotation changes) may differ on real spectra with near-ties.

## How it was run

```text
OPENBLAS_NUM_THREADS=1 simple_test_exec test=polar_matching_bench suite=<per_pair|replica|blocked|states_inner|matrix|matrix_t> nthr=11
OPENBLAS_NUM_THREADS=1 simple_test_exec test=polar_matching_bench suite=check nthr=11
```

Optional `nspace=`, `nstates=` (one case) and `nptcls=` (default 2000). Rows are appended to
`polar_matching_bench.tsv` and, for `check`, to `polar_matching_check.tsv` in the working directory. The program
is measure-only: no CTest entry, no thresholds. It is listed in section 1.3 of the test policy. Code:
`src/main/commanders/test/simple_polar_matching_bench.f90` (the layouts) and
`simple_commanders_test_highlevel_matchbench.f90` (cases, sweeps, checks, TSV).
