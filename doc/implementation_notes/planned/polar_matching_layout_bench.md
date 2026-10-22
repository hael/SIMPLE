# Polar matching layouts for many references: a measure-only benchmark

Planned 2026-10-21 for an unattended run on Hans's Mac (kit `~/agent_runs/runs/match_bench_mac`). It writes one
measure-only test program and a report; it changes no production code.

## Why

The dense probabilistic table (`refine=prob`) scores every sampled particle against every reference: `nspace`
projection directions times `nstates` states. On 2026-10-21 the maintainer moved the multi-state `prob_neigh` switch
back to stage 6, so the default five-stage multi-state `solve3D` runs the dense table in stages 3-5, stage 5 over
2,500 directions per state. The dense table is again the hot spot of multi-state runs. The question is how much
faster its scoring can be with a layout that reuses each reference across several particles, before anything is
changed in production.

## What a table entry costs today

For one particle and one reference, the euclid score over all in-plane rotations
(`gen_prob_likelihood_euclid_val` and `gen_corrs` in `src/main/pftc/simple_polarft_corr.f90`) is two sums over the
polar shells `k` of a particle-side spectrum times a reference-side spectrum: the particle's CTF² term against the
reference's power (`ft_ctf2 x ft_ref2`), and the CTF-weighted particle against the reference
(`ft_ptcl_ctf x conj(ft_ref)`). Each sum runs over `pftsz+1` angular frequency bins and is followed by one inverse
real FFT of length `nrots`. Since 2026-10-15 a particle's shift is applied on the particle side, so the references
stay fixed (`ft_ref`, `ft_ref2` are memoized per reference and half). The table fill (`fill_tab_range` in
`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`) gives each thread one particle at a time and
loops over all references.

An estimate at the sizes of Cyril's two-state set in stages 4-5 (crop 160, mask radius about 58 crop pixels, matching
limit about 6 Å: about 60 shells, `pftsz` about 182, `nrots` about 364): one pair reads about 176 KB of reference
spectra (two arrays) for about 180 kflop in the sums and about 15 kflop in the FFTs, about one flop per byte. With two
states and 2,500 directions the references take about 1.8 GB (both arrays, both halves), far more than any cache, so
every reference streams from memory again for every particle: the sums should be bound by memory bandwidth. The
two-state runs' stage 4 on the Dell (2,000 references, about 20,000 particles, about 115 s per iteration including
other work) works out at about 70 µs per pair and thread, close to streaming 176 KB at 2-3 GB/s per thread. This
benchmark replaces the estimate with measurements.

## The layouts compared

All four score the same synthetic data with the same in-plane sampling, at the same sizes:

1. **Per pair (today).** The production routine, called as `fill_tab_range` calls it: a thread takes a particle and
   scores it against every reference in turn.
2. **Blocked loop.** A thread takes a block of `B` particles and a tile of `R` references and, for each reference of
   the tile, scores every particle of the block with the same per-pair computation, summing over `k` in the same
   order. Each pair's sums are the same operations in the same order, so its distances must be byte-identical to
   layout 1. Measure `B` in {4, 8, 16, 32}; choose `R` so that a tile fits the core's L2 cache and report it.
3. **States inner.** For one direction, every state's reference is scored against the particle while the particle's
   data is in cache, with one batched inverse FFT per direction (an FFTW plan over `nstates` transforms). Report
   whether its distances are byte-identical to layout 1.
4. **Matrix product per angular bin.** For a tile of `R` references and a block of `B` particles, and for each angular
   frequency bin `j`, both sums become complex matrix products over the shells: references (R x nk) times particles
   (nk x B), through the BLAS `cgemm` SIMPLE links. This needs the references in a bin-major copy (`[j][k][r]`): build
   it once, report its size and build time. Then a batched inverse FFT per pair. Same sums in a different order: report
   the largest relative difference of the distances against layout 1, and the fraction of pairs whose best rotation
   differs.

The comparison of distances is made on the per-rotation distance vectors before the in-plane draw: the draw uses
random numbers whose order may change with the loop order, and it is the same code downstream in every layout.

## Sizes and cases

- Take the polar sizes from SIMPLE itself (`magic_pftsz` for the mask radius, the shell range for the limit) for box
  180 at 2.055 Å, crop 160, mask diameter 270 Å, limit 6 Å; record `pftsz`, `nrots` and the shell range used.
- Directions: 1,000 and 2,500. States: 1, 2 and 4. Particles: 2,000 per measurement (enough for stable timing),
  alternating between the two halves as production does. Objective euclid (production's for `solve3D`), with
  synthetic noise powers (sigma2) of realistic shape, CTF on, in-plane sampling as production's default
  (`prob_athres`).
- Synthetic spectra: random complex values of unit scale with a smooth radial falloff; no zeros or denormals.
  Timing does not depend on the values, but the arrays must be fully written before timing.
- Threads: 1, 5 (the performance cores) and 11 (all cores), dynamic schedule as production. Each timing is the best
  of three repetitions after one warm-up, with nothing else running on the machine (record the load).

## What is reported

For each layout, case and thread count: wall time, µs per pair per thread, pairs per second, the reference bytes
touched per second (an effective bandwidth), and the speedup against layout 1; for layouts 2 and 3 the identity
check, for layout 4 the differences. Also the memory each layout needs. A table in the report and a TSV file.

What the benchmark cannot say: the share of a real run's table time that scoring takes. The shift refinement of the
best references (L-BFGS from a coarse grid) depends on real data. The report names the measurement that settles it
(the table's timing profile in BENCH output of one real two-state run: seeding, scoring and refinement separately).

## The program

A measure-only program of `simple_test_exec` built like `cont_refine3D_bench` (commit 54c742363): no CTest entry, no
floors, listed in section 1.3 of `doc/policies/test_environment_policy.md`; `scripts/check_test_registry.py` must
still pass. The prototype kernels of layouts 2-4 live in the benchmark's own module on the test side and call the
production polar calculator for its setup and its memoized arrays; layout 1 calls the production routine itself. An
explicit interface to the BLAS `cgemm` goes in the benchmark module (no production file changes). Large objects are
`allocatable` (compile-time policy). Follow `simple_polarft_corr_tester.f90` for building a polar calculator from
parameters alone.

## Phases

### Phase 0: the harness and today's cost

1. Base Release build of the worktree (settings matched to `~/src/SIMPLE/build/CMakeCache.txt`; record the BLAS
   library CMake found, for example Accelerate or OpenBLAS) and the Debug build with the fast gate.
2. The program with layout 1 only: setup at the sizes above, timing loop, TSV output.
3. Layout 1 timed for every case and thread count.

Exit: layout 1's numbers for all cases, the sizes recorded, the fast gate passing.

### Phase 1: the other layouts

1. Layout 2, its identity check against layout 1 (all cases, at 1 and 11 threads, repeated five times at 11 threads),
   then its timings over `B`.
2. Layout 3, its identity check, its timings.
3. Layout 4, its differences, its timings over `R` and `B`.
4. Any parallel region added: tested at 11 threads, repeated, against one thread, with no string temporaries inside
   the region and work arrays sized by the team in use.

Exit: all four layouts measured and checked, or the ones that could not be made to work described with the reason.

### Phase 2: the report

The report `doc/implementation_notes/planned/polar_matching_layout_bench_report.md` in the worktree: the answer
first in plain words (how much faster scoring can be, with which layout, for one, two and four states, and at what
cost in memory and identity), then the table, then what it would take in production (the files and the in-plane
draw order question) and the measurement on real data that should come first. Then `full.diff` in the run
directory.

## The clock

The run has eight hours from launch. Phase 0 should take about two hours, Phase 1 about four, Phase 2 one. If time
runs short, layout 3 goes first, then the 5-thread measurements, then the 1,000-direction cases; the report is always
written.

## Common rules

- Nothing is committed; the run works in its own worktree; Hans's checkouts are read only.
- No production code changes. A defect found in production code is described in the report, not fixed.
- No floor or tolerance; this is a measurement. Identity checks are reported as they come out.
- Fortran rules: `use ..., only:`; THROW_HARD and THROW_WARN never continued with `//&`; no goto; no `.and.`
  short-circuit on allocatable arrays; no integer-overflow tricks; comments short.
