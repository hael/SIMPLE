# Production polar kernels: faster code for the same per-pair sums

Planned 2026-10-21 for an unattended night run on Hans's Mac (kit `~/agent_runs/runs/kernel_mac`). It follows the
measure-only benchmark of the same evening (`polar_matching_layout_bench.md` and its report, in the same directory
of the worktree) and works in the same worktree, where the benchmark program is present but not committed.

## Why

The benchmark scored synthetic polar spectra at the sizes of Cyril's two-state set (crop 160, 61 shells, `pftsz`
182, `nrots` 364). Two findings concern production code as it is, with no change of method:

1. **Production's per-pair sums are slower than the same source compiled elsewhere.** The "replica" (production's
   euclid sums of `gen_euclid_crvec`, unshifted branch, and its inverse FFT, compiled in the benchmark module on a
   memo the benchmark rebuilds) ran 1.37-1.43 times faster at one thread, 1.33-1.41 at five and 1.22-1.26 at eleven.
   In the Debug build the two give the same bits; in Release they differ in the last bit in about 4 % of the
   20-lowest lists, which the benchmark traced to how the compiler vectorizes the loop in each compilation unit.
   The suspected cause: the hot loop reads and writes arrays through components of the `polarft_calc` object
   (`self%ft_ref`, `self%ft_ctf2`, `self%crvec1(ithr)%c`, the last a pointer component shared with FFTW), so the
   compiler must assume they may overlap, and the accumulation mixes single and double precision per element
   (`cmplx(acc + wk*c)` with a double-precision `wk`).
2. **The distance conversion and the draw are a sixth of a pair.** At one thread a pair costs about 10.8 µs in the
   sums, 0.6 µs in the inverse FFT and 2.3 µs in turning 364 values into distances
   (`euclid_dist_from_crvec`, one function call per rotation, in double precision, with a branch) and drawing among
   the 20 lowest.

The dense table runs these kernels for every sampled particle against every reference. Since 2026-10-21 the
multi-state stages 3-5 use the dense table again, so a faster per-pair kernel speeds up every such stage, and every
single-state `prob` stage too, with no extra memory.

## Rulings that hold

From the maintainer, 2026-10-21: the speed changes stay; NU filtering at stage 5 stays; no change of resolution
policy. For this run: no change of method or layout (the matrix-product route waits for a measurement on the Dell);
production code changes only in the polar calculator.

## Items

1. **Diagnose.** Compile `src/main/pftc/simple_polarft_corr.f90` and the benchmark module with gfortran's
   vectorization reports (`-fopt-info-vec-optimized`, `-fopt-info-vec-missed`, in a scratch build only) and, if
   needed, read the generated assembly of the two inner loops. Record which loops vectorize in which unit and why
   not. This decides the shape of item 2.
2. **The euclid sums.** For `gen_euclid_crvec` (both branches) and `gen_euclid_crvec_bank`, move the accumulation
   over shells and angular bins into a private procedure of the submodule whose arguments are plain arrays
   (explicit shape or `contiguous`), with a local accumulator, called with the memo's sections. Keep the operation
   order and the precision of every element exactly as today (the same expression, the same conversions), so that
   the only difference can be the compiler's code generation. Then measure.
3. **The cc sums.** `gen_corrs` and `gen_corrs_bank` (the cc objective, used where `objfun=cc`) have the same
   structure. Extend the benchmark with a cc case if that is quick (test side only); otherwise apply item 2's pattern,
   check identity and correctness, and report it as not timed.
4. **The distance conversion.** In `sample_euclid_crvec`, `gen_prob_euclid_val`, `gen_best_euclid_val` and the
   bank path, compute the 364 distances in one loop that the compiler can vectorize (the same double-precision
   expression and the same treatment of overflow as `euclid_dist_from_crvec`, without a call per element). The
   result must be byte-identical. Leave the draw itself as it is unless an identical restructuring is obvious.

## How each item is judged

- **Speed.** The benchmark's `per_pair` layout (production) against its value on the base build, on all six cases
  of the benchmark, at 1, 5 and 11 threads, best of three after a warm-up, nothing else running. The target for
  item 2 is the replica's speed.
- **Identity.** The benchmark's check suite: the per-rotation distance vectors and the 20-lowest lists of production
  on the new build against the base build, at one thread, and each build at 11 threads repeated five times against
  its own one-thread run. Report the fraction of pairs that are byte-identical and the largest difference in units
  in the last place. Item 4 must be byte-identical.
- **Keep rule.** A change of items 2 and 3 is kept when it is at least 5 % faster at 11 threads on the largest case
  and its distances differ from the base build by at most two units in the last place; otherwise it is restored and
  the report says why. Byte identity is preferred: if a variant reaches the replica's speed with byte identity,
  keep that one.
- **Correctness.** The Debug build with the fast gate, the polar calculator's unit tests (`simple_polarft_corr_tester`
  and every unit suite that runs the polar calculator), and one high-level CTest entry that builds a probabilistic
  table (for example `simulated_workflow_1jxy`) passing its checks on the final build. Record whether that entry's
  outputs are byte-identical to the base build at one thread; a difference is not a failure if items 2-3 changed
  last bits, but it is reported.

## Phases

### Phase 0: base and diagnosis

Base Release build of the worktree as it is (the benchmark included), the benchmark's `per_pair` and `replica`
timings and the check suite on it (reuse the evidence of the benchmark run where it was made on the same tree and
machine; re-time `per_pair` once on the largest case to check drift), the Debug build with the fast gate. Item 1.

Exit: the base numbers, the diagnosis written in the journal.

### Phase 1: the kernels

Items 2, 3 and 4 in that order, each as its own item with a snapshot, a per-item diff, its measurement and its
identity check. Any parallel region touched: stress-tested at 11 threads against one thread.

Exit: every item kept or restored with its evidence; the correctness checks passing on the final build.

### Phase 2: the report

`doc/implementation_notes/planned/polar_kernel_codegen_report.md` in the worktree: the answer first (how much faster
the per-pair kernel is, whether it is byte-identical, what changed in the code), a table per item, the diagnosis
of why production was slower, and what it means for a real run (the share of scoring in a real table is still to
be measured on the Dell). Then the diffs: `full.diff` (everything since HEAD, the benchmark included) and
`kernel_only.diff` (this run's changes alone, against the snapshot taken before Phase 0).

### Phases 3 and 4: the x86 follow-up (Dell, run `kernel_x86`, 2026-10-22)

Phases 0-2 ran on the Mac and were committed as ee5b68be1. On the Dell (Intel Core Ultra 9 285, hybrid: cpu0-7
performance cores with their own L2, cpu8-23 efficiency cores in clusters of four sharing an L2; gfortran 15.2,
`-O3 -march=native`) the run `ms_prob5check` measured the new euclid kernel at 0.80-0.85x of the old one at 12 and
24 threads, 0.92x on one performance core and 0.81x on one efficiency core; cc 1.18x on a performance core, 0.92x
on an efficiency core. Its reports showed the old euclid loops already vectorized there with 32-byte (AVX2)
vectors and a run-time overlap check. Goal: understand why the new form is slower on x86 and find a form of the
sums that is at least as fast as the old one on the Dell, likely to keep the Mac's gain, and equal to the old
results to the last bit or two.

**Phase 3: diagnosis and candidates.** Builds of the old form (f3ccb5fbd's `simple_polarft_corr.f90` in an export
of HEAD) and HEAD. Diagnosis of the euclid and cc inner loops in both forms: gfortran's vectorization reports and
the generated assembly (vector width, shuffles for the complex interleaving, packing calls, the overlap check,
unrolling), `perf stat` if available. Candidates, each a complete form of the euclid sums (both branches of
`gen_euclid_crvec`, and `gen_euclid_crvec_bank`) and of the cc sums (`gen_corrs`, `gen_corrs_bank`), with every
element's expression and precision kept: the old inline loops; HEAD's explicit-shape helpers; the old loops with a
local `pointer, contiguous` alias of the accumulator (and of the memo slices if useful); helpers with assumed-shape
`contiguous` dummies; helpers accumulating into a plain per-thread work array copied into the FFT buffer once per
pair; any further form the diagnosis points to. Every candidate measured for euclid and cc on the two cases
1000x1 and 2500x4, pinned to cpu3 (performance) and cpu20 (efficiency), and at 12 and 24 threads, and its identity
against the old form. Exit: the diagnosis, the table of candidates and the identity results in the journal and in
Progress.

**Phase 4: the chosen form.** The chosen form in `simple_polarft_corr.f90`, built and checked, and the report
`polar_kernel_codegen_x86_report.md` (answer first; diagnosis; every candidate in every setting; identity; an
appendix with the full source of every candidate's sums and call sites, so the Mac can measure the same forms; a
commit message).

**Judging, fixed before any measurement.** A candidate is acceptable when, on the Dell, it is not slower than the
old form (f3ccb5fbd) by more than 2 % in any of the four settings (one performance core, one efficiency core,
12 threads, 24 threads) for euclid, and likewise for cc, and its distances differ from the old form's by at most two
units in the last place. Among acceptable forms prefer the one that tells the compiler that the accumulator and the
memo slices are contiguous and do not overlap (what the Mac needed). If no form is acceptable, the old form is
restored for that objective. Production changes only in `simple_polarft_corr.f90` and only to the form of the sums
and their callers' argument passing; the distance conversion of ee5b68be1 stays. Speed is judged on the Dell only:
it is closer to the machines SIMPLE normally runs on, and the Mac only has to build and pass its checks (maintainer,
2026-10-23).

### File table

| file | why | phase |
|---|---|---|
| `doc/implementation_notes/planned/polar_kernel_codegen.md` | this plan: Phases 3-4, file table, progress | 3 4 |
| `src/main/pftc/simple_polarft_corr.f90` | the chosen form of the euclid and cc sums | 4 |
| `src/main/commanders/test/simple_polar_matching_bench.f90` | benchmark, test side, only if the measurements need it | 3 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | benchmark driver, test side, only if the measurements need it | 3 4 |
| `doc/implementation_notes/planned/polar_kernel_codegen_x86_report.md` | the x86 report | 4 |

## The clock

The run has eight hours. Phase 0 about two, Phase 1 about four and a half, Phase 2 one. If time runs short, item 3
goes first; the report is always written.

## Common rules

- Nothing is committed; the run works in the worktree; Hans's checkouts are read only.
- Production changes only in the polar calculator (`simple_polarft_corr.f90`, and `simple_polarft_calc.f90` if an
  interface is needed). The benchmark may be extended on the test side.
- No floor or tolerance is relaxed; the keep rule above is fixed before the measurements.
- Fortran rules: `use ..., only:`; THROW_HARD and THROW_WARN never continued with `//&`; no goto; no `.and.`
  short-circuit on allocatable arrays; no integer-overflow tricks; comments short; dead code found in touched code is
  deleted.

## Progress

### Phases 3-4, 2026-10-23: the kernels on x86

#### Phase 3: diagnosis and candidates (Dell, 2026-10-22/23)

Done. Files changed: this plan (Phases 3-4, the file table, this section) and
`src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` (the benchmark's middle team is now
half the threads, so `nthr=24` times 24, 12 and 1 on the Dell; `nthr=11` still times 11, 5 and 1 on the Mac).
No production file changed in this phase. Journal: `~/agent_runs/kernel_x86/state/journal_phase_3.md`.

**Why the new form is slower on x86.** gcc 15.2 applies unroll-and-jam (factor 2) to the shell loop of the new
`euclid_sums` and to both shell loops of `cc_sums`. With plain array dummies it can prove that two shells may be
interleaved; through the old pointer accumulator it could not. On the Mac this gave two shells per pass with
de-interleaving loads and stores, a gain. On x86 with AVX2 gcc vectorizes the jammed loop with cross-lane permutes:
the euclid loop takes 276 instructions per 8 bins x 2 shells (17 per bin and shell; the old loop 15), among them
36 `vpermd`, 16 `vblendps`, 18 reloads of permute constants and 24 register spills to the stack; the cc loop has
20 `vpermd` per iteration. The clean loop the compiler also emits (12.5 instructions per bin and shell) runs only
for the odd last shell. Lane permutes are slow on the efficiency cores, hence 0.81x there against 0.92x on a
performance core. Neither `!$omp simd` on the bin loop nor `!GCC$ unroll 1` on the shell loop stops the jam
(checked in the compiler reports; not timed). The old loop was already vectorized at full AVX2 width; it loads and
stores the accumulator 8 bytes at a time, which costs less here than the permutes. A local `pointer, contiguous`
alias does not change gfortran's code: the accumulator stays strided. The distance conversion of ee5b68be1 is not
involved: it is byte-identical and alone makes the old sums 2.7 % faster on a performance core, 0.8 % on an
efficiency core. Writing two shells per pass by hand avoids the jam and its permutes and is the fastest form on
one thread, but it reads eight memory streams per pass (four arrays, two shells) and loses at 12 threads; forms
that read at most four streams per pass (one euclid shell, or two cc shells) keep their gain in a team.

**How it was measured.** Every candidate is HEAD's `simple_polarft_corr.f90` with only the sums and their callers
changed, built in one scratch tree (only that object recompiles), timed with the benchmark's production layouts
`per_pair` (euclid) and `per_pair_cc` on the cases 1000x1 and 2500x4. P1 and E1: one thread pinned to cpu3 and
cpu20, best of three after a warm-up, median of the pinned rounds (reproducible within 1 %). T12 and T24: the
unpinned teams of `nthr=24`, as production runs them; these vary 10-30 % between launches (placement of the team
on the two core types), so six rounds in rotating order (HEAD two, c14 four), median. Nothing else ran during the
timings; the machine's load was the benchmark itself.

Speed against the old form (old time / candidate time, above 1 is faster), 1000x1 / 2500x4; the old row gives
µs per pair and thread; team cells give the number of rounds. Candidates not timed in teams were already outside
the rule on one thread (c4, c5), share a loop with a timed one (c6 with c12/c13), equal the reference (c1), or
are diagnostics (d1); c7 is the old-loop version of the two-shell form that loses in teams.

euclid:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 14.39 / 14.42 | 23.64 / 23.65 | 24.10 / 30.34 (6) | 31.02 / 38.00 (6) |
| c1 old loops, new conversion | 1.028 / 1.029 | 1.008 / 1.008 | - | - |
| c2 HEAD helpers | 0.924 / 0.927 | 0.811 / 0.809 | 0.833 / 0.834 (2) | 1.019 / 0.829 (2) |
| c3 old loops, contiguous pointer alias | 1.023 / 1.025 | 1.002 / 1.001 | 0.964 / 0.969 (6) | 1.022 / 1.007 (6) |
| c4 helpers, contiguous assumed shape | 0.927 / 0.927 | 0.802 / 0.806 | - | - |
| c5 helpers, local work array | 0.929 / 0.929 | 0.814 / 0.801 | - | - |
| c6 helpers, two shells per pass | 1.155 / 1.150 | 1.092 / 1.096 | - | - |
| c7 old loops, two shells per pass | 1.090 / 1.077 | 1.063 / 1.055 | - | - |
| c9 helper call per shell | 1.123 / 1.127 | 1.119 / 1.120 | 0.997 / 1.042 (6) | 1.165 / 1.079 (6) |
| c12 euclid two shells, cc per shell | 1.154 / 1.151 | 1.103 / 1.101 | 0.892 / 0.959 (6) | 0.896 / 1.039 (6) |
| c13 euclid two shells, cc two shells per call | 1.156 / 1.153 | 1.039 / 1.033 | 0.897 / 0.924 (6) | 0.882 / 0.989 (6) |
| c14 euclid per shell, cc two shells per call | 1.148 / 1.150 | 1.127 / 1.126 | 0.999 / 1.062 (4) | 1.284 / 1.094 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 0.911 / 0.910 | 0.803 / 0.803 | - | - |

cc:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 17.61 / 17.68 | 20.39 / 20.50 | 25.51 / 34.25 (6) | 28.37 / 39.10 (6) |
| c1 old loops, new conversion | 0.998 / 1.001 | 1.002 / 1.001 | - | - |
| c2 HEAD helpers | 1.176 / 1.176 | 0.919 / 0.922 | 0.989 / 1.119 (2) | 0.984 / 1.110 (2) |
| c3 old loops, contiguous pointer alias | 1.182 / 1.182 | 1.282 / 1.281 | 0.947 / 1.135 (6) | 1.181 / 0.996 (6) |
| c4 helpers, contiguous assumed shape | 1.149 / 1.172 | 0.916 / 0.918 | - | - |
| c5 helpers, local work array | 1.183 / 1.183 | 0.908 / 0.907 | - | - |
| c6 helpers, two shells per pass | 1.410 / 1.427 | 0.997 / 0.996 | - | - |
| c7 old loops, two shells per pass | 1.055 / 1.059 | 1.013 / 1.010 | - | - |
| c9 helper call per shell | 1.350 / 1.341 | 1.383 / 1.386 | 1.110 / 1.172 (6) | 1.141 / 1.202 (6) |
| c12 euclid two shells, cc per shell | 1.347 / 1.339 | 1.370 / 1.380 | 1.263 / 1.146 (6) | 1.076 / 1.209 (6) |
| c13 euclid two shells, cc two shells per call | 1.386 / 1.363 | 1.297 / 1.284 | 1.472 / 1.254 (6) | 1.258 / 1.304 (6) |
| c14 euclid per shell, cc two shells per call | 1.387 / 1.398 | 1.305 / 1.307 | 1.491 / 1.218 (4) | 1.345 / 1.327 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 1.186 / 1.187 | 0.919 / 0.922 | - | - |

The tightest cell, euclid 1000x1 at 12 threads, was timed eight more times, interleaved, for old, c9 and c14:
c14 1.040x (24 threads 1.052x), c9 1.082x (1.133x); pooled with the six rounds c14 1.043x / 1.057x, c9 1.055x /
1.144x. The cc gain of every form that passes the memoized FT(REF) directly (c2, c3 and all helper forms) partly
comes from no longer copying it into `cmat2_many` when there is no shift.

**Identity against the old form** (the benchmark's `suite=identity`: 32 particles, the production kernels'
per-rotation vectors at one thread; for cc the decision of the Mac run holds: the judged distances are those the
table scores, 1 - cc, and the raw correlations are reported in ulp of each vector's largest |cc|, because they
cancel to about 1e-2). 2500x4 (320,000 pairs); 1000x1 gave the same bounds or tighter and no rotation change:

| form | euclid max ulp (dists / exp values) | euclid vectors identical | 20-lowest lists reordered | cc 1-cc max ulp | raw cc, ulp of vector max | cc vectors identical | 20-lowest lists reordered |
|---|---|---|---|---|---|---|---|
| c1 | 0 / 0 | 100 % | 0 | 0 | 0 | 100 % | 0 |
| c2 HEAD, c4, c5, d1 | 2 / 2 | 74.4 % | 1 | 1 | 4 | 77.3 % | 4 |
| c3 | 0 / 0 | 100 % | 0 | 1 (shifted only) | 6.5 | 100 % unshifted | 0 |
| c6 | 2 / 2 | 77.5 % | 4 | 1 | 4 | 77.5 % | 4 |
| c7 | 1 / 2 | 67.8 % | 6 | 1 | 2 | 100 % | 0 |
| c9, c14 | 2 / 2 | 67.8 % | 6 | 1 | 2 | 99.9 % | 0 |
| c12, c13 | 2 / 2 | 77.5 % | 4 | 1 | 2 | 99.9 % | 0 |

Two ulp occur only in the shifted-reference branch and in the bank's exp(-loss) values; the unshifted and bank
distances differ by at most one ulp. The best value of a pair never changes its rotation in any candidate. c14's
euclid bits equal c9's, its cc bits c13's.

**Reading.** Within 2 % of the old form in all four settings for both objectives and within two ulp: c9 and c14.
c3 misses (euclid 0.964 at 12 threads, cc 0.947); c12 and c13 miss (two-shell euclid, 0.88-0.96 in teams); HEAD
misses. Of the two, c14 is faster for cc in every setting (1.22-1.49x against 1.11-1.39x) and at least as fast for
euclid; both tell the compiler that the accumulator and the memo slices are contiguous and do not overlap
(explicit-shape helper dummies). Recommendation for Phase 4: c14, the euclid sums as one helper call per shell
(`euclid_shell_sum`) and the cc sums as one helper call per two shells (`cc_shells_sum`). c14 keeps unit-stride
access with no overlap check (what the Mac needed) and two shells per pass for cc; its euclid gives up the jam's two
shells per pass. The Mac's speed is not checked: speed is judged on the Dell.

**Exit items.** Diagnosis: above and in the journal, sections "Diagnosis" (compiler reports and assembly in
`scratch/diag/`, `scratch/diag/loopinfo/`). Table of candidates: above; raw logs `scratch/time/<round>/`, table
scripts `scratch/time/table.py` and `mdtable.py`, output `scratch/logs/mdtable.txt`. Identity: above; full
comparisons `scratch/logs/ident_table_1000x1.txt` and `_2500x4.txt`, comparer `scratch/keep/idcmp.c`; the old form's
reference vectors kept in `scratch/keep/ident_old/`. Candidate sources: `scratch/keep/cands/<form>/` (generator
`gen.py`), all under `~/agent_runs/kernel_x86/`.

#### Phase 4: the chosen form (Dell, 2026-10-23)

Done. Files changed: `src/main/pftc/simple_polarft_corr.f90` (the c14 form), this plan (this entry, the Progress
heading), the report `polar_kernel_codegen_x86_report.md`. Journal: `~/agent_runs/kernel_x86/state/journal_phase_4.md`.

- **Choice, by the rule, per objective.** Euclid: acceptable c9 and c14, which share the form (one helper call
  per shell, `euclid_shell_sum`). cc: acceptable c9 and c14; c14's form (one helper call per two shells,
  `cc_shells_sum`) is faster in every setting. Both forms pass plain array dummies, so the compiler knows the
  accumulator and the memo slices are contiguous and do not overlap. The combination is c14; no objective falls
  back to the old form. The file equals Phase 3's c14 except two doc-comment lines saying why the shell loop is
  split into helper calls.
- **Checks on the final tree** (Release build via `compile_clean.sh`, log `scratch/logs/build_repo_p4.log`): fast
  gate 15/15; `unit_pftc_align2D3D` 9/9 groups, 718/718 assertions; `unit_cart_align3D` 5/5 groups, 801/801
  assertions (`scratch/p4/*.log`); CTest `simulated_workflow_1jxy` passed, 17/17 metrics (`scratch/p4/ctest_1jxy.log`,
  `metrics_1jxy.tsv`); benchmark `suite=check` at 24 threads: every layout thread-reproducible, production at 24
  threads equal to one thread (`scratch/p4/check.log`, `check_table.md`); `suite=identity` of the repository build
  against the old form gives exactly c14's bits, and its 24-thread repeats equal one thread
  (`scratch/p4/ident_final_vs_old.txt`).
- **Re-timed against the old form** (one pinned round, two team rounds in opposite order, machine otherwise idle;
  `scratch/logs/retime_table.txt`), speed 1000x1 / 2500x4: euclid P1 1.150 / 1.148, E1 1.126 / 1.125, T12 1.129 /
  1.050, T24 1.245 / 1.074; cc P1 1.376 / 1.395, E1 1.323 / 1.308, T12 1.113 / 1.117, T24 1.170 / 1.288. All cells
  above 1, confirming Phase 3.
- **Report** written with the answer, the diagnosis, every candidate in every setting, identity, the Phase 4
  checks, decisions, deviations, open points, a commit message and the full source of every form.

Exit items: choice (above); checks (fast gate, both unit suites, the simulated workflow, `suite=check`, re-timing:
all passed); report and this Progress entry written. Paths are under `~/agent_runs/kernel_x86/`.
