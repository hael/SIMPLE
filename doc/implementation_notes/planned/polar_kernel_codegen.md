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
