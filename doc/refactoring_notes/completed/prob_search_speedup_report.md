# The probabilistic search speed-up: report

Report of the interactive session of 14 October 2026 that carried out
`prob_search_speedup.md` (beside this file). The code was written in a session
that could not compile or run anything; the maintainer built it, ran the gates
and the cases on the bgal benchmark data, and reported back. Nothing is
committed. The full record, phase by phase, is in section 11 of the plan.

## Terms

- **The pre-step**: the probabilistic table that `refine=prob` and
  `refine=prob_neigh` build before every refine3D iteration. A driver
  (`prob_align`, `prob_align_neigh`) launches parts (`prob_tab`,
  `prob_tab_neigh`); each part scores its particles against the references and
  refines the shifts of the best ones; the driver reads the parts' tables and
  assigns every particle a state, a projection and an in-plane angle.
- **Score loop**: scoring one particle against every reference over all
  in-plane rotations. **Refinement**: the shift search (a 5 x 5 coarse grid,
  then L-BFGS-B) on the best references of a particle. **Seed**: the shift
  found once per particle on its previous best reference, used to score every
  reference in the single-state table.
- **Reference frame / particle frame**: a shift applied to the reference in
  its own polar frame, or to the particle. They differ by the in-plane
  rotation; a particle-frame shift is the same at every rotation, so the
  particle can be shifted once and scored against unshifted references.
- **Thread-seconds**: time summed over threads, as opposed to wall time.

## What changed

**Benchmark files (exact).** `PROB_TAB_BENCH_ITERnnn_PARTp.txt` and
`PROB_ALIGN_BENCH_ITERnnn.txt` time the pre-step: wall-clock steps, a
thread-seconds block (seed, score loop, refinement, load balance) and the
search objects' counters. The neighbourhood table (`PROB_TAB_NEIGH_BENCH`) now
writes the same thread-seconds block and the number of references scored
(written, not yet run).

**`prob_fast` (approximate, developer key, default off; refine3D, solve3D,
solve3D_chunks, solve3D_addon).** One key switches on
three changes, each in its own scope, so the same key serves single-state and
multi-state runs:
- *particle-frame seed*, single-state table: the seed is applied to the
  particle once and every reference is scored from that one shifted copy;
- *particle-frame grid*, multi-state table: the 5 x 5 coarse shift grid is a
  lattice of shifted particles, built once per particle and shared by all its
  refined references;
- *quadratic step*, multi-state table: where the grid minimum has a full 3 x 3
  neighbourhood and the fit is well conditioned and beats the grid minimum, a
  quadratic fit replaces the L-BFGS-B search; otherwise L-BFGS-B runs as
  before. With the key on, the multi-state fill uses dynamic OpenMP
  scheduling, because the particles' costs then differ.

The new kernels (`memoize_shifted_ptcl`, `gen_objfun_vals_bank`,
`gen_likelihood_val_bank`) include the angular Nyquist bin, so a score is
exact against the true shifted loss; today's shifted-reference kernels drop
that bin and differ slightly (maintainer's ruling). The three changes started
as three keys and were merged into `prob_fast` at the end (ruling 7).

**Subspace map of the neighbourhood search (exact, default path).**
`build_general_tbox` built the map from subspace to full projection
directions by calling `sym_dists` for every (subspace, full) pair. The
subspace count grows with `nspace`, so the build grew with its square: 3.5 s
per process at `nspace` 5000 (d2), paid by the `prob_align_neigh` driver, each
`prob_tab_neigh` part and each refine3D part. `sym%sym_dist_table` makes the
symmetry copies once per direction and the table from dot products, with the
same operands, branches and tie handling (ruling 8; written, not yet run).

**A crash fixed (default path).** A state that received no sampled particles
on iteration 1 of a run started from input volumes stopped the run with
`fsc_state02.bin does not exist`: volassemble dropped the state, but the
refine3D driver still counted it as populated (the unsampled particles keep
their labels) and read its FSC. The driver now keeps the input volume as that
state's reference and skips its FSC, plot, registration and postprocess
(written, not run; the maintainer waived the test).

## Evidence

Single state, `solve3D`, bgal (4551 particles), `prob_seed_ptcl=yes` (now
`prob_fast=yes`) against keys off:

| | keys off | key on |
| --- | --- | --- |
| whole run | 631.1 s | 505.9 s |
| FSC 0.5 / 0.143 | 4.29 / 3.75 A | 4.29 / 3.75 A |
| `prob` pre-step, iterations 1-8 | 198.1 s | 82.0 s |
| score loop, thread-s | 3462 | 907 |

Two states, `refine3D`, 2 iterations, `update_frac=0.2`, `lp=8`, `vol1`/`vol2`
the half-maps of the single-state volume (5000 references, 910 particles per
iteration):

| | keys off | grid + fit | grid + fit, dynamic schedule |
| --- | --- | --- | --- |
| table step, iterations 1 / 2 | 34.4 / 33.5 s | 16.3 / 15.3 s | 15.2 / 14.4 s |
| refinement, thread-s (part 1) | 141 | 40 | 42 |
| gradient evaluations per search | 12.2 | 6.1 | 6.1 |
| load balance | 0.98 | 0.84-0.97 | 0.98-0.99 |
| whole run | 96.9 s | 60.4 s | 58.9 s |
| score, iterations 1 / 2 | 0.637 / 0.548 | 0.644 / 0.556 | 0.649 / 0.551 |

State populations agreed within sampling (814/96 against 815/95 and 808/102 in
iteration 1). The small second state's resolution varied between runs (13.1 to
29.7 A at FSC 0.143) with a cFAR near zero; a repeat of the keys-off run, to
measure that spread, was not made.

Where the single-state run's 505.8 s go (key on): `prob` pre-step 82 s,
`prob_neigh` pre-step 160 s, refine3D parts and merge 143 s, strategy setup
39 s, in-iteration assembly 26 s, final reconstruction 48 s, start and end
7 s. The subspace-map fix targets the `prob_neigh` pre-step and the refine3D
parts of the `nspace` 5000 iterations; expected about 10 s per such iteration.

## Decisions

| Decision | Where |
| --- | --- |
| New kernels include the angular Nyquist bin | plan, section 11 (ruling of 14 October) |
| One key, `prob_fast`, for the three approximations | plan, ruling 7 |
| Subspace-map table on the default path | plan, ruling 8 |
| Blocked table loop (Phase 4) skipped | plan, section 11 |
| Edge-centred stencil for the quadratic step rejected | plan, section 11 |

**Rejected changes.**
- *Phase 4, the blocked table loop.* Its own criterion was not met (score loop
  28 thread-s against refinement 40 with the keys on), and the loop is
  compute-bound: 1.93e-5 s per particle-reference pair with 6 threads against
  2.44e-5 s with 24 threads on 20 cores, a difference the size of the
  oversubscription. Blocking could recover about 1 s of a 12 s fill. Not
  written.
- *A 3 x 3 stencil shifted inside the lattice for an edge grid minimum.* 41% of
  the multi-state searches fell back because the grid minimum lay on the edge
  of the lattice. The shifted stencil converted 14% of them; the rest failed on
  curvature or step instead, because an edge minimum mostly means the optimum
  lies beyond the lattice. No measurable gain, and 0.23 px of fit error in the
  unit gate. Reverted; the diff is `prob_search_speedup_edge_stencil_rejected.patch`
  beside this file.

## Tests

- `unit_pftc_align2D3D`: frame conversions, the shifted-particle bank against
  the direct sum (2e-7) and against today's kernel, the particle-frame grid and
  the quadratic step on the phantom. The CC reference of the bank test was
  corrected to the unweighted FFT objective; passes.
- `3D probability table` (new suite, `simple_eul_prob_tab_tester.f90`): the
  shift recorded for an assigned, unrefined candidate; passes.
- `symmetry` suite, `test_sym_dist_table` (new): the table against
  `sym_dists` pair by pair in nine point groups, and the same nearest
  neighbours; not yet run.

## Open

- **Not yet built or run:** the merge into `prob_fast`, `sym_dist_table` and
  its test, the neighbourhood-table instrumentation, and the crash fix. The
  next build needs the unit gate (`unit_ori`, `unit_pftc_align2D3D`) and one
  single-state `solve3D` with `prob_fast=yes`, compared with
  `bench/7_solve3D`: same FSC, `startup/setup` of the neighbourhood parts down
  from 3.5 s, iterations 12-17 shorter.
- **`prob_fast` is a developer key pending promotion.** Promotion makes the
  three changes the only path and removes the key (release 4 keeps no backward
  compatibility). A multi-state `solve3D` with the key on has not been run.
- **The particle-frame seed in the neighbourhood table.** In single-state
  `prob_neigh` every reference of a particle is scored with the same seed, the
  pattern the key exploits in the full table. Worth building only if the
  instrumented bench file shows a large score loop.
- **The spread of the small state's resolution** between identical runs.

Follow-up work belongs in a new planned note.

## Skill text

None proposed: nothing here changes how code is written in the repository.
