# Polar kernels on x86: the report

Report of Phases 3 and 4 of the plan `polar_kernel_codegen.md` (same directory), the x86 follow-up of the Mac run
committed as ee5b68be1. Unattended run `kernel_x86` on the Dell (Intel Core Ultra 9 285: 8 performance cores with
their own L2, 16 efficiency cores in clusters of four sharing an L2; gfortran 15.2.1, `-O3 -march=native`, system
OpenBLAS 0.3.15, FFTW 3.3.5), 2026-10-22/23. The only production file changed is
`src/main/pftc/simple_polarft_corr.f90`. Nothing is committed.

## The answer

- **Why x86 was slower.** With the sums in helpers with plain array arguments, gcc can prove that two shells may be
  processed in one pass and applies unroll-and-jam to the shell loop. On the Mac that gave two shells per pass with
  de-interleaving loads and stores, a gain. On x86 with AVX2 gcc vectorizes the jammed loop with lane permutes
  (36 `vpermd` per pass), reloads of permute constants and register spills: 17 instructions per bin and shell
  against 15 for the old loop. Permutes are slow on the efficiency cores most of all, hence HEAD's euclid kernel at
  0.81x of the old one there and 0.92x on a performance core. The old loops could not be jammed, because their
  accumulator is a pointer component; they were already vectorized at full AVX2 width.
- **The form kept (c14).** The euclid sums call a helper once per shell (`euclid_shell_sum`), the cc sums call one
  once per two shells (`cc_shells_sum`). The helpers still take plain array arguments, so the compiler knows unit
  stride and no overlap (what the Mac needed), but there is no shell loop around a bin loop in one procedure left
  to jam. Each element's expression and precision are unchanged; the distance conversion of ee5b68be1 stays.
- **Its speed against the old form (f3ccb5fbd) here:** euclid 1.15x on a performance core, 1.13x on an efficiency
  core, 1.00-1.28x in teams of 12 and 24 threads; cc 1.39x and 1.31x alone, 1.22-1.49x in teams. Against HEAD the
  euclid kernel is 1.24x faster on a performance core and 1.39x on an efficiency core.
- **Identity.** Distances differ from the old form's by at most one ulp, two in the shifted-reference branch and
  in the bank's exp(-loss) values; the best rotation of a pair never changes; 6 of 320,000 lists of the 20 lowest
  rotations change order (HEAD: 1 for euclid, 4 for cc). cc distances (1 - cc): at most one ulp, 99.9 % of the
  vectors byte-identical.
- **The Mac.** Speed is judged on the Dell, which is closer to the machines SIMPLE normally runs on; the Mac only
  has to build and pass its checks (maintainer, 2026-10-23). Whether c14 keeps the Mac's euclid gain of ee5b68be1
  (its euclid now runs one shell per pass there too) is not measured.

## What changed, phase by phase

- **Phase 3** (`review/phase_3.diff` under `~/agent_runs/kernel_x86`): the plan gained Phases 3-4, the judging rule,
  a file table and the Phase 3 progress (diagnosis, every candidate's timings, identity). The benchmark driver
  `simple_commanders_test_highlevel_matchbench.f90` times its middle team at half the threads (24 -> 24, 12, 1 on
  the Dell; 11 -> 11, 5, 1 on the Mac as before). No production change.
- **Phase 4** (`review/phase_4.diff`): `simple_polarft_corr.f90` gets the c14 form: `euclid_sums` calls the new
  `euclid_shell_sum` once per shell; `cc_sums` calls the new `cc_shells_sum` once per two shells (once more for an
  odd last shell). Callers, argument passing and the distance conversion are HEAD's. This report and the plan's
  Progress complete the phase.

## Why HEAD was slower here

Read from gfortran's vectorization and loop reports (`-fopt-info-vec-all`, `-fopt-info-loop-optimized`) and the
assembly of the hot loops, with the build's flags:

- The old euclid loops vectorize 4 complex values per iteration (32-byte vectors), behind a run-time overlap check.
  The accumulator `crvec1(ithr)%c` is a pointer component of unknown stride, so it is loaded and stored 8 bytes at
  a time: 61 instructions per 4 bins.
- HEAD's `euclid_sums` contains a clean loop of 50 instructions per 4 bins, but it runs only for the odd last
  shell. For more than two shells the code enters the unroll-and-jammed version (gcc reports "applying unroll and
  jam with factor 2" for the shell loops of `euclid_sums` and both of `cc_sums`): 276 instructions per 8 bins x
  2 shells, of them 36 `vpermd` (cross-lane permutes), 16 `vblendps`, 18 reloads of permute constants and 24 stack
  accesses for spilled registers. The cc loop has 20 `vpermd` per iteration; the old cc loop none.
- No packing copies (`_gfortran_internal_pack`) in either form's hot path; the FFTW buffer is 32-byte aligned, so
  HEAD's vector stores never split a cache line; the distance conversion vectorizes without masked stores and is not
  involved (byte-identical; alone it makes the old sums 2.7 % faster on a performance core).
- What does not help: `!$omp simd` on the bin loop and `!GCC$ unroll 1` on the shell loop leave the jam in place;
  a local `pointer, contiguous` alias of the accumulator does not change gfortran's code at all (still strided).
- Writing two shells per pass by hand removes the permutes and is the fastest form on one thread (c6, c12, c13:
  euclid 1.15x), but it reads eight memory streams per pass (four arrays, two shells) and falls to 0.88-0.96x at
  12 threads. The likely reason, not proven: cores that share an L2 run out of hardware prefetch streams. Forms
  that read at most four streams per pass (one euclid shell, or two cc shells of two arrays) keep their gain in a
  team. This is why c14 treats the two objectives differently.

## The candidates in every setting

How measured: each form is HEAD's file with only the sums and their callers replaced, built in one scratch tree;
the benchmark's production layouts `per_pair` (euclid) and `per_pair_cc` on the cases 1000x1 (1,000 references)
and 2500x4 (10,000). P1 and E1: one thread pinned to cpu3 (performance) and cpu20 (efficiency), best of three after
a warm-up, median of two or three pinned rounds (repeatable within 1 %). T12 and T24: unpinned teams as production
runs them; they move 10-30 % from launch to launch, so six rounds in rotating order (c14 four, HEAD two), median.
Nothing else ran during the timings. Speed = old time / form time, 1000x1 / 2500x4; the old row gives µs per pair
and thread; team cells give the number of rounds. Forms without team timings were already outside the rule on one
thread (c4, c5), share their euclid loop with a timed form (c6 with c12 and c13), equal the reference (c1), or are
diagnostics (d1); c7 is the old-loop version of the two-shell form.

euclid:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 14.39 / 14.42 | 23.64 / 23.65 | 24.10 / 30.34 (6) | 31.02 / 38.00 (6) |
| c1 old loops, new conversion | 1.028 / 1.029 | 1.008 / 1.008 | - / - | - / - |
| c2 HEAD helpers | 0.924 / 0.927 | 0.811 / 0.809 | 0.833 / 0.834 (2) | 1.019 / 0.829 (2) |
| c3 old loops, contiguous pointer alias | 1.023 / 1.025 | 1.002 / 1.001 | 0.964 / 0.969 (6) | 1.022 / 1.007 (6) |
| c4 helpers, contiguous assumed shape | 0.927 / 0.927 | 0.802 / 0.806 | - / - | - / - |
| c5 helpers, local work array | 0.929 / 0.929 | 0.814 / 0.801 | - / - | - / - |
| c6 helpers, two shells per pass | 1.155 / 1.150 | 1.092 / 1.096 | - / - | - / - |
| c7 old loops, two shells per pass | 1.090 / 1.077 | 1.063 / 1.055 | - / - | - / - |
| c9 helper call per shell | 1.123 / 1.127 | 1.119 / 1.120 | 0.997 / 1.042 (6) | 1.165 / 1.079 (6) |
| c12 euclid two shells, cc per shell | 1.154 / 1.151 | 1.103 / 1.101 | 0.892 / 0.959 (6) | 0.896 / 1.039 (6) |
| c13 euclid two shells, cc two shells per call | 1.156 / 1.153 | 1.039 / 1.033 | 0.897 / 0.924 (6) | 0.882 / 0.989 (6) |
| c14 euclid per shell, cc two shells per call | 1.148 / 1.150 | 1.127 / 1.126 | 0.999 / 1.062 (4) | 1.284 / 1.094 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 0.911 / 0.910 | 0.803 / 0.803 | - / - | - / - |

cc:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 17.61 / 17.68 | 20.39 / 20.50 | 25.51 / 34.25 (6) | 28.37 / 39.10 (6) |
| c1 old loops, new conversion | 0.998 / 1.001 | 1.002 / 1.001 | - / - | - / - |
| c2 HEAD helpers | 1.176 / 1.176 | 0.919 / 0.922 | 0.989 / 1.119 (2) | 0.984 / 1.110 (2) |
| c3 old loops, contiguous pointer alias | 1.182 / 1.182 | 1.282 / 1.281 | 0.947 / 1.135 (6) | 1.181 / 0.996 (6) |
| c4 helpers, contiguous assumed shape | 1.149 / 1.172 | 0.916 / 0.918 | - / - | - / - |
| c5 helpers, local work array | 1.183 / 1.183 | 0.908 / 0.907 | - / - | - / - |
| c6 helpers, two shells per pass | 1.410 / 1.427 | 0.997 / 0.996 | - / - | - / - |
| c7 old loops, two shells per pass | 1.055 / 1.059 | 1.013 / 1.010 | - / - | - / - |
| c9 helper call per shell | 1.350 / 1.341 | 1.383 / 1.386 | 1.110 / 1.172 (6) | 1.141 / 1.202 (6) |
| c12 euclid two shells, cc per shell | 1.347 / 1.339 | 1.370 / 1.380 | 1.263 / 1.146 (6) | 1.076 / 1.209 (6) |
| c13 euclid two shells, cc two shells per call | 1.386 / 1.363 | 1.297 / 1.284 | 1.472 / 1.254 (6) | 1.258 / 1.304 (6) |
| c14 euclid per shell, cc two shells per call | 1.387 / 1.398 | 1.305 / 1.307 | 1.491 / 1.218 (4) | 1.345 / 1.327 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 1.186 / 1.187 | 0.919 / 0.922 | - / - | - / - |

Compiled but not timed (their code is HEAD's or c6's, checked in the assembly): c8 (c6 with cc as explicit bin
loops: cc still jammed), c10 (`!$omp simd`), c11 (`!GCC$ unroll 1`).

The tightest cell for the kept form, euclid 1000x1 at 12 threads, was timed eight more times interleaved with old
and c9: c14 1.040x (24 threads 1.052x); pooled with the six rounds 1.043x (1.057x). The cc gain of all forms that
pass the memoized FT(REF) directly (c2, c3 and every helper form) partly comes from no longer copying it into
`cmat2_many` when there is no shift.

**The rule.** Acceptable (within 2 % of the old form in all four settings, both objectives, and within two ulp):
c9 and c14. c3 misses (euclid 0.964 and cc 0.947 at 12 threads), c12 and c13 miss (two-shell euclid in teams),
HEAD misses. Chosen per objective: euclid, one helper call per shell (c9 and c14 share it); cc, two shells per
helper call (c14's, faster than c9's in every setting). Both tell the compiler that the accumulator and the memo
slices are contiguous and do not overlap. Together that is c14.

## Identity

The benchmark's `suite=identity`: 32 particles, the production kernels' per-rotation vectors at one thread,
compared bit for bit with the old form's. For cc the judged distances are those the table scores, 1 - cc (the Mac
run's reading); the raw correlations cancel to about 1e-2 and are reported in ulp of each vector's largest |cc|.
2500x4 (320,000 pairs); 1000x1 gave the same bounds or tighter and no reordered list.

| form | euclid max ulp (dists / exp values) | euclid vectors identical | 20-lowest lists reordered | cc 1-cc max ulp | raw cc, ulp of vector max | cc vectors identical | 20-lowest lists reordered |
|---|---|---|---|---|---|---|---|
| c1 | 0 / 0 | 100 % | 0 | 0 | 0 | 100 % | 0 |
| c2 HEAD, c4, c5, d1 | 2 / 2 | 74.4 % | 1 | 1 | 4 | 77.3 % | 4 |
| c3 | 0 / 0 | 100 % | 0 | 1 (shifted only) | 6.5 | 100 % unshifted | 0 |
| c6 | 2 / 2 | 77.5 % | 4 | 1 | 4 | 77.5 % | 4 |
| c7 | 1 / 2 | 67.8 % | 6 | 1 | 2 | 100 % | 0 |
| c9, c14 | 2 / 2 | 67.8 % | 6 | 1 | 2 | 99.9 % | 0 |
| c12, c13 | 2 / 2 | 77.5 % | 4 | 1 | 2 | 99.9 % | 0 |

Two ulp occur only in the shifted-reference branch and in the bank's exp(-loss) values. The best value of a pair
never changes its rotation in any form. c14's euclid bits equal c9's and its cc bits c13's.

## Checks of the final tree (Phase 4)

On a Release build of the repository with the c14 form (`compile_clean.sh`, `-O3 -march=native`):

| check | result |
|---|---|
| fast gate | 15/15 passed |
| `simple_test_exec test=unit_pftc_align2D3D` | 9/9 groups, 718/718 assertions |
| `simple_test_exec test=unit_cart_align3D` | 5/5 groups, 801/801 assertions |
| CTest `simulated_workflow_1jxy` | passed (785 s), 17/17 metrics pass |
| benchmark `suite=check`, 24 threads | every layout's 24-thread repeats equal their own one-thread run; production at 24 threads equals one thread |
| benchmark `suite=identity` against the old form | exactly c14's bits (table above); five 24-thread repeats of all four kernels equal one thread, both objectives and cases |
| re-timing against the old form | below |

Re-timing of the repository build against the old form (one pinned round, two team rounds in opposite order,
nothing else running), speed = old time / new time, 1000x1 / 2500x4:

| objective | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| euclid | 1.150 / 1.148 | 1.126 / 1.125 | 1.129 / 1.050 | 1.245 / 1.074 |
| cc | 1.376 / 1.395 | 1.323 / 1.308 | 1.113 / 1.117 | 1.170 / 1.288 |

Evidence, all under `~/agent_runs/kernel_x86/scratch/`: `logs/build_repo_p4.log`, `p4/*.log`, `p4/check_table.md`,
`p4/metrics_1jxy.tsv`, `p4/ident_final_vs_old.txt`, `logs/retime_table.txt`; Phase 3's timings `time/<round>/`,
tables `logs/mdtable.txt`, identity `logs/ident_table_*.txt`, compiler reports and assembly `diag/`.

## Decisions taken under delegation

All recorded in the journals (`~/agent_runs/kernel_x86/state/journal_phase_3.md`, `journal_phase_4.md`) and, where
they shape the result, in the plan's Progress.

- **Waiting for an idle machine.** Another run's jobs held all 24 cores until 02:35; no timing ran before 02:47.
  Builds, compiler reports and the identity runs (bits, not time) ran before, niced.
- **The benchmark's middle team** is half the threads (test side), so that `nthr=24` times the 12-thread team the
  run's instructions ask for; the Mac's `nthr=11` still gives 11, 5, 1.
- **cc identity** judged on the distances the table scores, 1 - cc, as the Mac run did; raw correlations reported
  in ulp of each vector's largest |cc|.
- **Team timings** judged on the median of six unpinned rounds (they move 10-30 % between launches), plus eight
  interleaved launches for the closest cell. Production runs teams unpinned, so pinning them was not an option.
- **Candidates beyond the run's list.** The diagnosis pointed at the jam, so forms that avoid it were added (c6-c14);
  c8, c10 and c11 were dropped after their compiler reports showed HEAD's or c6's code.
- **Which forms got team timings:** those not already ruled out on one thread or equal to a timed form (see the
  table's note).
- **The choice** made per objective as the rule says; the two winners combine into the complete form c14.
- **A Release build** of the repository (`compile_clean.sh`, fast gate included) rather than the debug script, so
  that the re-timing compares like with like; the unit suites and the simulated workflow ran on it.
- **Disk:** 48 GB of identity vectors deleted after their comparisons were written to text; the old form's
  vectors kept under `scratch/keep/ident_old`.

## Deviations from the plan

- The plan (the Mac run's) had no Phases 3-4 and no file table; both were added from the run's instructions, with
  the judging rule fixed before any measurement.
- The benchmark driver change (test side) is the only change outside the plan's production file and documents.
- The re-timing of the chosen form is one pinned round and two team rounds against the old form, not one run, for
  the same noise reason as above.
- `simulated_workflow_1jxy` was checked for passing (17/17 metrics), not compared byte for byte with a run of the
  old build; the identity suite covers the kernels themselves.

## Open for the maintainer

- **Real runs.** All numbers are from the synthetic benchmark. What share of a real dense table's time the kernels
  take on the Dell was measured by `ms_prob5check` for HEAD (score loop 52-68 % of the fill time); a real run with
  c14 was not made.
- **The memory-stream explanation** of the two-shell euclid's team loss is a hypothesis: `perf` is not installed on
  the Dell, so there are no counter measurements.
- **Not touched:** `gen_raw_euclid_vals_impl` still accumulates through the pointer component (as the Mac report
  already noted); it is not on the dense table's path.

## Commit message

```text
pftc: per-shell helper calls for the euclid and cc sums

With the shell loop and the bin loop in one helper, gcc unroll-and-jams two
shells; on x86 (AVX2) the jammed loop needs lane permutes and register spills,
and the euclid kernel ran at 0.81-0.92x of the old inline loops. The euclid
sums now call a helper once per shell and the cc sums once per two shells, still
with plain array arguments (unit stride, no overlap). On the Dell euclid is
1.13-1.15x and cc 1.31-1.39x the old form on one core, and at least as fast in
teams; distances change by at most 2 ulp. The benchmark's middle team is half
the threads.
```

## Appendix: the source of every form

The full [candidate source archive](supporting/polar_kernel_codegen_x86_candidates.md)
preserves the original appendix. The choice, timings, identity checks and limitations remain above.
