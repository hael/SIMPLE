# Production polar kernels: the report

Report of the plan `polar_kernel_codegen.md` (same directory). Unattended run on Hans's Mac (Apple M3 Pro,
5 performance and 6 efficiency cores, 18 GB), 2026-10-21/22. The only production file changed is
`src/main/pftc/simple_polarft_corr.f90`. Nothing is committed.

## The answer

- **The per-pair euclid kernel is 1.30 times faster at 11 threads, 1.34 at 5 and 1.36 at 1**, in all six cases of
  the benchmark. On the largest case (2,500 directions, 4 states) a pair costs 23.8 µs per thread instead of 31.0
  at 11 threads. That is the speed of the "replica" (22.9 µs), the target the plan set.
- **The cc kernel is 1.50 times faster at 11 threads, 1.61 at 5 and 1.70 at 1.**
- **The learning: production was slow because of one pointer.** The sums accumulate into the per-thread FFT
  buffer `self%crvec1(ithr)%c`. It is a pointer component, so the compiler cannot know that its elements are
  adjacent. It therefore stored the results one complex number at a time and shuffled lanes to do so. Passing the
  same arrays to a small procedure with plain array arguments tells the compiler that the data are contiguous and
  do not overlap. The arithmetic in the source is unchanged.
- **It is not byte-identical, only last-bit identical.** The euclid distances differ from the old build by at most
  one unit in the last place (ulp) without a shift and with the particle bank, and two ulp with a shifted
  reference. The best value of a pair is unchanged in 99.8 % of pairs, its rotation in all of them. The cc
  distances the table uses (1 − cc) differ by at most one ulp. The reason is code generation, not arithmetic: in
  the Debug build, where nothing is vectorized, the new kernel gives exactly the replica's bits, which were already
  known to equal production's there. No byte-identical variant was found.
- **The distance conversion is now vectorized and byte-identical.** It gains only 0.2-1.3 %: in that part of a
  pair, the selection of the 20 lowest distances and the draw cost more than the conversion.
- **A real run is unchanged at the level the tests see.** `simulated_workflow_1jxy` passes on the old and the new
  build with the same 17 metrics, and all 206 maps and stacks it writes are byte-identical.

## Terms

- **Pair**: one particle scored against one reference over all 364 in-plane rotations: the sums over 61 shells
  and 182 angular bins, one inverse FFT, the conversion to 364 distances, the draw among the 20 lowest.
- **µs/pair/thread**: wall time × threads / pairs. Lower is better.
- **Speedup**: old time / new time, same case and threads.
- **ulp**: unit in the last place of a single-precision number; "1 ulp" is the smallest possible difference.
- **Base**: the worktree before this run, with the benchmark's new test-side suites (production code unchanged).
- **Replica**: production's euclid sums compiled inside the benchmark module on a rebuilt copy of the memo (from the
  previous run). It was 1.35 times faster than production, which is what started this work.
- **Byte-identical**: every bit of every value the same.

## What changed in the code

All in `simple_polarft_corr.f90`, the polar calculator's objective submodule.

| item | change |
|---|---|
| 2, euclid sums | new private `euclid_sums(np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc)` with explicit-shape arguments, called by `gen_euclid_crvec` (both branches) and `gen_euclid_crvec_bank` with the memo's sections and `crvec1(ithr)%c`. The expression and the precision of every element are the old ones. |
| 3, cc sums | new private `cc_sums(n1, ldb, nk, a, b, l_conj, acc)`, called by `gen_corrs` and `gen_corrs_bank`. Without a shift, `gen_corrs` now passes the memoized `FT(REF)` directly instead of copying it into the `cmat2_many` buffer first. |
| 4, distances | new private `euclid_dists_from_crvec(n, r, norm, dists)` replaces the per-element function `euclid_dist_from_crvec` (deleted) in `gen_best_euclid_val`, `gen_prob_euclid_val`, `sample_euclid_crvec` and `gen_prob_power_euclid_val`. Same double-precision expression and threshold; every element is converted (the value clamped first, so a huge value is never converted to single precision) and those above the threshold are then set to `huge`, which lets the loop vectorize. |

No interface, type component, parallel region or caller outside the submodule changed. The benchmark (test side)
gained `suite=identity` (the production kernels' full distance vectors, written to a file per build and compared
between builds and between team sizes) and `suite=per_pair_cc` (the cc objective's timing).

## Why production was slower (item 1)

From gfortran's vectorization reports and the assembly of both compilation units (Release flags, scratch build):

- Production's euclid loop was vectorized two complex numbers per iteration. Every result was loaded and stored
  with separate 8-byte accesses at a run-time stride, because `crvec1(ithr)%c` is a pointer. The memo arrays were
  loaded as vectors and rearranged with lane shuffles. An aliasing check ran before every shell.
- The replica's loop, the same source with the accumulator a dummy argument, used de-interleaving loads and stores
  (`ld2`/`st2`), four complex numbers per iteration, and handled two shells per pass over the accumulator.
- The mixed precision (single-precision products, a double-precision weight, rounding to single on every
  accumulation) is the same in both, so it was not the cause.
- The distance conversion was not vectorized at all: its `if`/`else` is control flow the vectorizer does not
  handle. A `merge()` did not help either (gfortran keeps it as a branch); "convert, then overwrite" did.
- The cc sums: one complex number per iteration in `gen_corrs`' first sum, the cross term not vectorized.

## Item 2: the euclid sums

Speed, `per_pair` (production's `gen_likelihood_val`, as the dense table calls it), best of three after a warm-up,
2,000 particles (200 at one thread). µs/pair/thread, base → item 2 (speedup):

| directions | states | 11 threads | 5 threads | 1 thread |
|---|---|---|---|---|
| 1000 | 1 | 30.97 → 23.93 (1.29x) | 20.16 → 14.98 (1.35x) | 18.70 → 13.90 (1.35x) |
| 1000 | 2 | 31.01 → 23.98 (1.29x) | 20.04 → 14.98 (1.34x) | 18.73 → 13.90 (1.35x) |
| 1000 | 4 | 30.94 → 24.01 (1.29x) | 19.96 → 15.04 (1.33x) | 18.86 → 13.91 (1.36x) |
| 2500 | 1 | 31.02 → 24.01 (1.29x) | 19.89 → 15.08 (1.32x) | 18.78 → 13.94 (1.35x) |
| 2500 | 2 | 31.04 → 24.05 (1.29x) | 19.99 → 15.00 (1.33x) | 18.95 → 13.92 (1.36x) |
| 2500 | 4 | 30.97 → 24.05 (1.29x) | 20.07 → 15.12 (1.33x) | 18.97 → 14.02 (1.35x) |

The replica on the largest case: 22.86 / 14.70 / 13.73. The keep rule (at least 5 % at 11 threads on the largest
case, at most two ulp) is met: 29 %.

Identity against base. 64 particles × 2,000 references (1000x2) and 32 × 10,000 (2500x4), one thread:

| kernel | pairs byte-identical | largest difference |
|---|---|---|
| unshifted (`gen_euclid_crvec`, the dense table's path) | 63 % | 1 ulp (1.2e-7 absolute) |
| shifted reference | 0.2 % | 2 ulp (1000x2), 1 ulp (2500x4) |
| particle bank (`gen_euclid_crvec_bank`), distances | 3 % | 1 ulp |
| particle bank, `exp(-loss)` values | 0.7 % | 2 ulp |
| 20 lowest distances and their rotations | 96 % of lists | 1 ulp; rotation order differs in 2 of 80,000 lists |
| best value (`gen_best_euclid_val`) | 99.8 % | 1 ulp; the best rotation never changes |

One variant was tried to reach byte identity: assumed-shape `contiguous` arguments instead of explicit shape. It
gave exactly the same bits, so the compiler generates the same code. The per-element operations in the new
assembly look the same as the old; the difference is how the 182 bins are split between the vector loop and its
remainder. The source cannot steer that while keeping the same expression.

## Item 3: the cc sums

Speed, `per_pair_cc` (`gen_likelihood_val` with `objfun=cc`), timed on the smallest and the largest case:

| directions | states | 11 threads | 5 threads | 1 thread |
|---|---|---|---|---|
| 1000 | 1 | 28.79 → 19.75 (1.46x) | 18.80 → 11.65 (1.61x) | 17.54 → 10.35 (1.69x) |
| 2500 | 4 | 30.37 → 20.70 (1.47x) | 19.99 → 12.34 (1.62x) | 17.78 → 10.44 (1.70x) |

Identity against base:

- The cc distances the table scores (1 − cc, clamped to [0, 1]) differ by at most 1 ulp; 8 % of the vectors are
  byte-identical, 69 % of the 20-lowest lists.
- The raw correlations (`gen_objfun_vals`, `gen_objfun_vals_bank`) differ by at most 1.1e-8 absolute, which is
  5 ulp of each vector's largest |cc|. On these synthetic data |cc| < 0.045: the sums cancel, so a last-bit change
  of a sum is large relative to the small result. The rotation of the highest correlation never changes.
- Decision recorded during the run: the keep rule's "distances" were taken as the distances the table scores
  (1 − cc), which meet the two-ulp limit. If Hans means the raw correlations, item 3 does not meet the rule as
  measured here and should be reverted (its diff stands alone: `review/item_03_ccsums.diff`).

## Item 4: the distance conversion

The first version (two `merge()` calls) was byte-identical but gained nothing (0.99-1.01x): the reports showed
the loop still not vectorized. The kept version converts every element and then overwrites those above the
threshold. It vectorizes (16-byte vectors, unroll 4) and is byte-identical to item 3 in every kernel, objective and
case. A standalone check against the old function on 2,000,003 values, including the threshold's neighbourhood,
±huge, ±infinity, NaN and ±tiny, found no differing bit.

Speed against item 2: 0.999-1.013x, median about 1.006x. In the replica's profile (previous run) the conversion and
the draw together were 2.3 µs of a 13.8 µs pair at one thread; the vectorized conversion removes only part of the
conversion's share. The selection of the 20 lowest (a heap) and the draw now dominate that part.

## All items together

Final build against base, `per_pair`, every case: 1.297-1.307x at 11 threads, 1.331-1.353x at 5, 1.354-1.367x at 1.
`per_pair_cc` largest case: 1.50x / 1.61x / 1.70x. Identity against base: as item 3 for both objectives.

## Correctness

| check | result |
|---|---|
| every build, every kernel, 11 threads × 5 against one thread | byte-identical (both objectives, both cases) |
| Debug build (bounds-checked) and fast gate, base and final | 15/15 passed both times |
| Debug build, `suite=identity` and `suite=check` at 11 threads | clean; the replica equals production exactly |
| Release `unit_pftc_align2D3D` (polar calculator, in-plane, probability tables) | 9/9 groups, 698/698 assertions |
| Release `unit_cart_align3D` | 5/5 groups, 788/788 assertions |
| `simulated_workflow_1jxy` (CTest, Release, one thread), base and final | passed; 17 metrics identical; all 206 `.mrc`/`.mrcs` byte-identical |

In the simulated workflow 31 of 292 output files differ: timings (BENCH files), absolute paths (the two builds sit
in directories of different names; the project files differ already at extraction, before any polar code runs),
and one 8-byte header field of `sigma2_state.bin`, consistent with its layout digest. The sigma2 values are equal.
The last-bit differences of the kernels did not change any map of this workflow.

No parallel region was added or changed; the stress tests above cover the kernels the regions call.

## What it means for a real run

The dense probability table runs these kernels for every sampled particle against every reference, so a stage
that builds a dense table scores pairs about 1.3 times faster on this Mac. What share of a real table's time is
scoring is still to be measured: the table's timing profile (`ts_score` against `ts_refine` in
`prob_fill_profile`) of one real two-state `solve3D` run on the Dell. The Dell may also differ: the previous run
estimated that it is limited by memory bandwidth, not arithmetic, in which case faster arithmetic gains less there.

## Not done, and observations

- Not changed, outside this run's scope: `gen_raw_euclid_vals_impl` (single-precision weights) has the same loop
  shape over `crvec1(ithr)%c` and would gain from the same treatment. The shift-matrix products in the shifted
  branches were not looked at.
- Dated documents (`doc/implementation_notes/completed/drop_legacy_box_division.md`) and the generated
  `doc/code_overview/fortran-indexes/*` still name `euclid_dist_from_crvec`; the indexes regenerate.
- No defect was found in the code touched.

## Caveats

- One machine, one compiler: gfortran 16.2 (`-O3 -mtune=generic`, OpenMP), FFTW 3.3.11, OpenBLAS 0.3.27,
  `OPENBLAS_NUM_THREADS=1`. Another compiler or version may vectorize differently, in either direction, and the
  last bits will differ again.
- Synthetic data. Speed does not depend on the values. The identity statistics (how often a best rotation or a
  20-lowest list changes) may differ on real spectra with near-ties.
- Background during the timings: an editor and the agent (up to about a third of a core), backup and security
  services (a few %). Base timings matched the previous run's within 1 %.

## Evidence and how to reproduce

Run directory `~/agent_runs/runs/kernel_mac`: `journal.md` (every step and decision), `evidence/` (TSVs, tables,
logs), `review/item_NN_*.diff` (one diff per item), `full.diff` and `kernel_only.diff`.

```text
OPENBLAS_NUM_THREADS=1 simple_test_exec test=polar_matching_bench suite=per_pair nthr=11
OPENBLAS_NUM_THREADS=1 simple_test_exec test=polar_matching_bench suite=per_pair_cc nspace=2500 nstates=4 nthr=11
OPENBLAS_NUM_THREADS=1 simple_test_exec test=polar_matching_bench suite=identity nspace=1000 nstates=2 nptcls=32 nthr=11
```

`suite=identity` writes `polar_kernel_identity_<euclid|cc>_<case>.bin` (the one-thread vectors) and
`polar_kernel_identity.tsv` (five runs at the full team against one thread). Two builds' files were compared with a
small script that counts byte-identical pairs and the largest difference in ulp.

## Commit message

```text
pftc: vectorize the per-pair polar sums and the distance conversion

Move the euclid and cc sums of the per-pair kernels into private procedures
with plain array arguments. The accumulator is a pointer component, which kept
the compiler from using unit-stride vector code. The per-pair euclid kernel is
about 1.3x faster, cc about 1.5x, with the same expressions; distances change in
the last bit (at most 2 ulp). The euclid distance conversion is now one
vectorized loop, byte-identical. The polar matching benchmark gains an identity
suite that compares the production kernels between builds, and a cc timing.
```
