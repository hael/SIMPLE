# Global sigma2 only: review report

Report of the unattended run that removed `sigma_est=group` from SIMPLE,
2026-10-12 to 2026-10-13, for review of the uncommitted diff. The plan, with
the full record of every phase in its section 13, is beside this file:
`global_sigma2_only.md`. The per-phase diffs are in
`/home/elmlundho/agent_runs/global_sigma2_only/review/phase_N.diff` (N = 0 to
4). Logs, benchmark files and run scripts are under
`/home/elmlundho/agent_runs/global_sigma2_only/scratch/`.

## Terms

- **sigma2**: the noise power spectrum, one value per Fourier shell
  (resolution ring). The Euclidean objective (`objfun=euclid`) and the
  maximum-likelihood regularised reconstruction (`ml_reg=yes`) weight every
  shell by its inverse.
- **Half**: the even or the odd half-set of the particles. Each half is
  refined and reconstructed on its own, so the two maps can be compared.
- **Record**: one particle's residual spectrum, the power of the difference
  between the particle image and its matched reference projection.
- **The store**: `sigma2_state.bin`, the file that holds every particle's
  record and the two global spectra (one per half), each the mean of the
  valid records of that half's active particles. It is updated as a
  transaction: the workers write range files, the master builds a
  *candidate* file and *publishes* it by an atomic rename.
- **Rebuild (bootstrap)**: seeding a missing or unusable store from the
  particles' image power spectra (`calc_pspec`).
- **FSC**: Fourier shell correlation; resolutions are quoted where it falls to
  0.5 and to 0.143.
- **PCG**: preconditioned conjugate gradient, the iterative reconstruction
  backend (`rec_backend=pcg`).
- **FLEX**: SIMPLE's flexibility analysis (`flex_pca`, `cls_expansion`).
- **Workflows**: `solve3D` determines a 3D map de novo (from class averages
  or particles) and `refine3D` refines one; `solve2D` and `refine2D` are the
  2D classification counterparts; `make_cavgs` writes class averages;
  `reconstruct3D` reconstructs from fixed orientations. Each iteration of a
  matcher (the alignment search) updates the store.
- **Distributed / shared memory**: a run split into parts, each a separate
  process (`nparts`), with a master that merges their output; or one process
  with threads.
- **Generation**: the store's commit counter; every update publishes the
  next generation.

## What changed, and why

Before this run, `sigma_est=group` gave each particle stack its own noise
spectrum per half. Every workflow now uses one global spectrum per half,
which is what every workflow already used by default. The per-particle
records stay in the store, because a fractional update replaces only some
of them and the global spectra must reflect every particle's latest
residual. Release 4 keeps no backwards compatibility, except that `.simple`
project files of earlier releases still read; a store of the old kind is
rebuilt.

### Phase 0: timings and baseline (`review/phase_0.diff`)

Benchmark lines for the store transaction were added in SIMPLE's benchmark
format: the candidate preparation (refine2D and refine3D strategies), the
four steps of the commit (range merge, reduction, validation re-read,
publication), and the workers' record read and range write. The cases of
the plan were run on the unchanged code to give the baseline: beta-
galactosidase, 5,513 particles, box 256 at 1.275 Å per pixel; `solve3D_cavgs`
(S), `solve3D` (A), distributed and shared-memory `refine3D` (B, C), PCG
`reconstruct3D` (D), `solve2D` (E), `make_cavgs` (F), and a store written
under `sigma_est=group` (G) to test the rebuild.

### Phase 1: one rebuild entry, the parameter removed (`review/phase_1.diff`)

- Every workflow that reads or updates the store now goes through one
  rebuild entry (`ensure_sigma2_for_iteration`) and one identity check
  (`canonical_sigma2_consumable`). Before, five workflows had their own copy
  of the check and `refine2D` and `make_cavgs` had none.
- `sigma_est` and its derived flag are gone from the parameters, the user
  interface, every commander, the stream, the `solve3D` manifest (within
  manifest version 2) and the living documents.
- On copies of G, `refine2D` and `make_cavgs ml_reg=yes` refuse the old
  store, rebuild it and complete.

### Phase 2: the store and the in-memory model (`review/phase_2.diff`)

- The store accepts only one global spectrum per half. The per-stack
  grouping, the STAR-seed provenance and the direct range mode are refused,
  and the rebuild entry then rebuilds. The STAR converter `sigma2_convert`
  is retired.
- In memory, `euclid_sigma2` holds two spectra instead of a table with one
  spectrum per particle. Every consumer reads its particle's half: the
  polar and Cartesian matchers, gridding, class averaging, PCG (one spectrum
  per half and a scalar weight per particle), FLEX and `cls_expansion`.
- The renames of the plan's section 5, without aliases (`commit_sigma2`,
  `load_sigma2`, `read_global`, ...).
- Six defects of the existing code were fixed under the maintainer's
  rulings: shared-memory `refine3D` with `mkdir=yes` nesting its directory;
  a distributed master waiting forever on dead parts; a segmentation fault
  in the FLEX phantom test fixture; the `solve3D_addon` test refusing its
  own command line; the distributed PCG solve stopping on an empty state;
  and `cls_expansion` writing memberless sub-class halves as active blank
  averages.

### Phase 3: kernel and reconstruction optimisations, bootstrap (`review/phase_3.diff`)

Each optimisation was measured against the build without it, on the case
that exercises it, alternating the two builds.

- Kept: the shell weights k/sigma2(k) and each particle's weighted
  self-term computed once per spectra load; a cached reconstruction-plane
  grid with the per-half inverse spectra (`simple_fplane_grid.f90`); class
  bucketing in class-average restoration; the residual path without copies;
  a shell lookup table and reusable slots in the Cartesian calculator.
- The GPU path follows the per-half weights; there is no GPU on the test
  machine to measure it.
- Dropped: folding the weight into the particle transforms (not measurably
  faster; the patch is kept under `scratch/keep/bench/phase3/O7_folding/`).
- The bootstrap's workers write per-half sums instead of every row
  (`simple_sigma2_pspec_sums.f90`); the master writes the records and sums
  the global spectra while writing. The stored spectra are bit-identical,
  and the normalisation is preserved exactly, including the average image
  divided by the number of active particles in both halves.

### Phase 4: the streaming commit, documents and close (`review/phase_4.diff`)

- **The streaming commit.** Before, one matcher iteration copied the
  committed store into the candidate, merged each worker's range file into
  it (opening and syncing the whole file for each), re-read it to form the
  global spectra, and re-read it again to validate them. Now the strategy
  only validates the committed store and removes files a crashed
  transaction of the same generation left behind
  (`sigma2_state_prepare_update`). The commit
  (`sigma2_state_assemble_candidate`) creates a fresh candidate and fills it
  in blocks of 4,096 rows through one file handle: each row comes from the
  range file that covers it, or else from the committed store. The valid
  rows of active particles are summed per half as they pass, the global
  spectra are written, and the file is synced once before publication. The
  commit still refuses an uncovered or doubly covered row, a range of
  another transaction, an empty half and non-finite global spectra. The
  range files stay until publication, so a failed commit can be rerun.
- The workers take the transaction's identity (next generation, particle
  layout) from the committed store's header instead of from a prepared
  candidate.
- **The block-copying merges.** The chunk merge and `merge_projects`
  (`simple_projfile_utils.f90`) used to read each source store whole into
  memory, write it through its own open and sync, then reduce and validate
  by re-reading. They now call `sigma2_state_concatenate`, which copies in
  the same blocks, sums while copying and syncs once.
- **The validation re-read moved to the tests** (ruled). The reduction that
  re-reads every record (`sigma2_state_reduce`) stays as the reference the
  tests compare against and that test fixtures use to make a store. The
  old merge, validation and commit routines are deleted.
- **Tests** (`simple_sigma2_state_tester.f90`): the streamed global spectra
  equal a full reduction, on 4 and on 9,000 rows (three blocks); rows outside
  the schedule are taken from the committed store; an uncovered, a doubly
  covered row and a range of another generation are refused with the
  committed store unchanged byte for byte; a publication that fails leaves
  the committed store and the range files, and the commit then reruns from
  the same ranges; the concatenation keeps the row order and its global
  spectra equal a full reduction. The Cartesian calculator's tester follows
  the new commit.
- **Documents.** The grouped-sigma wording is gone from the `solve3D`,
  `solve3D_cavgs`, `refine3D` and PCG reconstruction policies; the stale
  grouping comments are gone from the code, and so is a `rec3D_backends`
  test description that still spoke of symlinked sigma2 group files.
  History documents stay as written; this report and the completed plan
  are the completion note.

## Evidence against Phase 0

All timed runs used the Phase 0 command lines on an out-of-tree Release
build, one at a time on the otherwise idle 24-core machine, with
`nparts=4 nthr=6` (distributed) or `nthr=24` (shared memory). Logs:
`scratch/logs/case_<name>.log`; benchmark files:
`scratch/keep/bench/<phase>/<case>/`.

### Accuracy, time and memory

Acceptance: every final FSC 0.143 resolution within the larger of the
Phase 0 spread and 3 % of its Phase 0 value; for `solve2D`, the maintainer's
rule (mean score at least 0.377, particle-weighted mean class resolution at
most 11.49 Å).

| Case | Phase 0 | Phase 3 (all optimisations) | Phase 4 (streaming commit) |
| --- | --- | --- | --- |
| A1 `solve3D` | 15:00, 3.19 GB, 4.35 / 3.84 Å | 15:00, 3.19 GB, 4.35 / 3.84 Å | 15:11, 3.21 GB, 4.35 / 3.89 Å |
| B1 `refine3D`, distributed | 3:10, 3.03 GB, 8.82 / 4.41 Å | 3:16, 3.03 GB, 8.82 / 4.35 Å | 2:58, 3.03 GB, 8.82 / 4.41 Å |
| C `refine3D`, shared memory | 2:11, 3.70 GB, 8.82 / 4.60 Å | 2:16, 3.69 GB, 8.82 / 4.60 Å | 2:04, 3.69 GB, 8.59 / 4.53 Å |
| D PCG `reconstruct3D` | 1:43, 6.59 GB, 3.80 Å | 1:48, 6.11 GB, 3.80 Å | not repeated (ruled) |
| E1 `solve2D` | 3:13, 0.41 GB, score 0.389, 10.85 Å | 3:07, 0.41 GB, 0.389, 10.96 Å | 2:49, 0.43 GB, 0.389, 11.19 Å |

Phase 0 spreads between two runs: A 0.06 / 0.00 Å, B 0.00 / 0.06 Å. A2, B2
and E2 were also run in Phases 2 and 3 (plan, section 13). F (`make_cavgs`)
completed in every phase up to 3 and was not repeated in Phase 4 (ruled).

### The store transaction

Seconds, summed over the iterations of the run (43 for A1, 5 for B1 and C,
27 for E1); for the workers, the slowest partition of each iteration. The
store has 5,513 rows of 128 shells (2.82 MB); a `refine3D` iteration
updates 1,103 rows (`update_frac=0.2`), a `solve2D` iteration all of them.

| Line | A1 P0 | A1 P4 | B1 P0 | B1 P4 | C P0 | C P4 | E1 P0 | E1 P4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| candidate preparation | 0.641 | 0.000 | 0.031 | 0.000 | 0.024 | 0.000 | 0.145 | 0.000 |
| merge + reduction + validation (P0) / streaming assembly (P4) | 0.846 | 0.251 | 0.104 | 0.029 | 0.049 | 0.028 | 0.572 | 0.155 |
| publication | 0.254 | 0.265 | 0.027 | 0.030 | 0.031 | 0.027 | 0.174 | 0.177 |
| commit total, inside the commander | 1.093 | 0.525 | 0.130 | 0.061 | 0.082 | 0.057 | 0.744 | 0.339 |
| commit, as the strategy sees it | 1.512 | 0.933 | 0.245 | 0.180 | 0.167 | 0.135 | 0.964 | 0.546 |
| workers' record read | 0.037 | 0.000 | 0.005 | 0.000 | 0.006 | 0.005 | 0.016 | 0.000 |
| workers' range write | 1.225 | 1.270 | 0.146 | 0.087 | 0.084 | 0.079 | 0.224 | 0.312 |

The workers' range write is unchanged code; its variation is the machine's.
On these few megabytes the whole transaction was already under 1 % of an
iteration in Phase 0 (about 0.3 % for `refine3D`, 0.7 % for `solve2D`); the
gain matters on large stores (see the sparse ranges below).

### Full passes over the store per matcher iteration

Counted from the code paths each iteration runs (store-sized reads and
writes; header and global-section reads are a few hundred bytes):

| | Reads | Writes | File syncs (nparts = 4) |
| --- | --- | --- | --- |
| Before | 5: workers' record read; candidate copy; range files in the merge; reduction re-read; validation re-read | 3: workers' range write; candidate copy; merged rows | 7: candidate copy, four range merges, global spectra, publication |
| After | 2: workers' record read; range files in the assembly | 2: workers' range write; candidate | 2: candidate, publication |

The workers' range files are each synced once, before and after. The merges
of chunk and project stores go from one full read and write per source
plus two re-reads of the result, with a sync per source, to one read and
one write with a single sync.

### Simulated workflow tests

Release build, `ctest -R '^<test>$' -V`. Phase 4 ran them side by side
after the timed cases (ruled).

| Test | Phase 0 | Phase 3 | Phase 4 |
| --- | --- | --- | --- |
| `simulated_workflow_1jxy` (truth correlation, FSC 0.143) | 0.899, 4.46 Å | 0.901, 4.54 Å | 0.902, 4.46 Å |
| `simulated_workflow_6vxx` (recorded, not judged) | 0.877, 5.43 Å | 6 of 8 runs pass | 0.880, 5.31 Å (1 run, passes) |
| `cont_refine3D_1jxy` (median rotation error, stages a/b/c) | 0.91 / 0.23 / 0.91° | 0.91 / 0.23 / 0.91° | 0.91 / 0.23 / 0.91° |
| `solve3D_addon` (union map against truth, FSC 0.143) | failed on the unchanged code | 0.938, 5.03 Å | 0.944, 4.93 Å; rerun 0.941, 4.93 Å |

`solve3D_addon` became testable in Phase 2 (one of the six defects fixed);
its baseline is Phase 2's 0.949 and 4.83 Å (floors 0.90 and about 6 Å).
Its closest check is the excess of the new cohort's pose error over the
frozen particles', limit 5°: 0.50° in Phase 2, 2.83° in Phase 3, and 4.76°
and 3.78° in two Phase 4 runs of the same build. The second run was made
because the first came close to the limit. The two runs differ by 1°, so the
measure varies from run to run.

### Gates

Every phase ended with `./compile_debug.sh` free of compiler warnings and
the fast gate passing 15 of 15 suites, `scripts/check_test_registry.py`
passing, a clean Release build and no file mode changes. Phase 4: the
same, plus the unit suites of the store, the Cartesian calculator, the
project and reconstruction (`scratch/logs/p4_debug_build_final.log`,
`p4_unit_final.log`, `p4_registry.log`).

## Decisions taken under delegation

Each is recorded in the plan's section 13 under its phase unless noted.

1. **Command lines of the baseline** (Phase 0): `ncls=100` for `solve2D`
   (the user guide's example); `calc_pspec` run distributed; case C with
   `mkdir=no` and F in a copy of E1's final directory, because the unchanged
   code failed both ways otherwise (two defects, fixed later under a
   ruling). Later phases repeat these lines unchanged.
2. **The final reconstruction keeps its own decision** (Phase 1; plan
   sections 3 and 7): an unusable store sends it to `bootstrap_rec3D`, which
   seeds and runs a residual pass at the final sampling. A plain rebuild
   would skip that pass.
3. **Where the rebuild check runs** (Phase 1): masters and shared-memory
   runs only, not the parts of a distributed reconstruction; new at the entry
   of `refine2D` and of `make_cavgs ml_reg=yes`. `refine2D` now writes its
   even/odd split before building references, so a `make_cavgs` inside its
   initialisation seeds on the same halves.
4. **A Phase 1 defect fixed in Phase 2**: the rebuild entry read the store
   path recorded in the project, so a new execution directory consumed the
   previous run's store. It now resolves the store against the working
   directory, as the consumers do.
5. **Half lookup** (Phase 2): a particle with `eo` outside 0 and 1 is refused
   (`sigma2_half`), also inside parallel loops; the polar path alone maps
   `eo = -1` to the even half, as before. Rebuilding a refused store restarts
   at generation 1.
6. **Reading of the first `solve2D` ruling** (Phase 2): the mean class
   resolution taken as the log's particle-weighted average. The maintainer
   confirmed this reading in the second ruling.
7. **Memberless sub-class halves** (Phase 2, third ruling): the work project
   of `refine3D_states_cavgs` had to learn the new per-half state flags,
   although the ruling expected no downstream change. Recorded in the
   `refine3D_states_cavgs` policy.
8. **Keep or drop each optimisation by measurement** (Phase 3): the table in
   section 13; folding dropped.
9. **The plane cache's key** (Phase 3): the padded box and the lowest shell,
   with the cache dropped whenever the spectra are reloaded, instead of the
   store generation the plan named (`euclid_sigma2` does not hold it, and
   FLEX's unit spectra have none).
10. **The Cartesian residual keeps its allocating interface** (Phase 3): it is
    formed from the calculator's own double-precision sums, twelve tester call
    sites rely on it, and the saved copy is negligible.
11. **The shape of the streaming commit** (Phase 4, this report and section
    13): the candidate is created by the commit, not before the workers;
    the strategy's preparation validates and clears stale files; the workers
    read the transaction identity from the committed header; the commit's
    benchmark lines become "streaming assembly", "publication" and "total".
12. **What stays of the reduction** (Phase 4): `sigma2_state_reduce` is kept
    as the full re-reading reduction that the tests compare the streamed
    spectra against and that fixtures use; `sigma2_state_validate_science`,
    `sigma2_state_merge_local_ranges`, `sigma2_state_commit`, the copying
    branch of `sigma2_state_create_candidate` and its compatibility check are
    deleted, having no other user.

## Deviations from the plan

- **Six defects of the existing code fixed** (Phase 2), by the maintainer's
  rulings; they were outside the plan. Three stops preceded them.
- **The `solve2D` acceptance replaced** by the maintainer's rule, after the
  unchanged code itself missed the plan's class-count and median rule.
- **Living documents' `sigma_est` lines moved to Phase 1** (from Phase 4),
  because Phase 1's exit required the tree to hold no `sigma_est`.
- **The plane cache key** and **the Cartesian residual interface**: decisions
  9 and 10 above. **The GPU weights** are coded but not measured.
- **Phase 4 validation trimmed** by ruling: A, B, C and E once each, D and F
  not repeated, `simulated_workflow_6vxx` recorded and not judged,
  `flex_pca_blobs` not rerun. One extra `solve3D_addon` run was made for
  information, because the first came close to a limit.
- **File-table additions** beyond the plan, each with its reason in the
  table: the `cls_expansion` commander (Phase 1); the queue controller, FLEX
  and addon test files of the rulings (Phase 2); the plane-cache and
  bootstrap-sums modules (Phase 3); the Cartesian calculator's tester, the
  `solve3D` manifest's comments and the `rec3D_backends` test description
  (Phase 4).
- **The skill file is not edited**: its proposed text is below.

## Proposed text for the skill file

In `.github/skills/simple-refine3d/references/bayesian-model.md`, section
"Regularization and Statistical Controls", replace the two bullets

```
- `objfun=euclid`
  Noise-normalized Euclidean objective with sigma estimation.
- `sigma_est`
  Controls sigma calculation policy.
```

with

```
- `objfun=euclid`
  Noise-normalized Euclidean objective. The noise model is one spectrum per
  even/odd half: the mean of the per-particle residual spectra in the
  project's committed sigma2 store (`sigma2_state.bin`). Each matcher
  iteration updates the records of the particles it evaluated and commits a
  new generation; a missing or incompatible store is rebuilt from particle
  power spectra at workflow entry (`ensure_sigma2_for_iteration`).
```

## Left open

- **The sparse ranges.** Workers would write only the rows they updated, and
  the streaming commit would take the others from the committed store, so a
  fractional update with fraction f would read the store about once and write
  it about 1 + f times per iteration instead of twice each. The plan makes it
  depend on the maintainer's ribosome measurement (751,729 rows, on a spinning
  disk): if, after the streaming commit, the transaction is still a visible
  share of an iteration. The lines to read are the strategy files' "sigma2
  candidate prep" and "sigma2 commit", the commit file's "streaming assembly"
  and "publication", and the matcher files' "sigma2 record read" and "range
  write".
- **`simulated_workflow_6vxx`, accepted as flaky for now.** Its runs are not
  reproducible (multithreaded reduction order changes the 2D step from its
  first iteration). Counts: Phase 0 1 of 1 pass; Phase 2 2 of 2; before
  Phase 3, 7 of 8; after Phase 3, 6 of 8; Phase 4, 1 of 1. Every
  failing run lands on the same wrong solution: whole-volume correlation 0.51
  to 0.55 against the truth and a truth FSC 0.143 of 49.9 to 83.2 Å, where
  passing runs reach 0.87 to 0.89 and 4.6 to 6.1 Å. Settling it needs about 20
  runs per build and, if the rates differ, bisecting the Phase 3 changes that
  alter numbers. By Phase 3's measurements only the reconstruction-plane cache
  changes results (the last bits of gridding and restoration); the shell
  weights, class bucketing, PCG and the bootstrap are bit-identical.
- **`flex_pca_blobs`, a known failure of the existing code.** All its FLEX
  checks pass (label accuracy 0.9995 and 0.9975, state maps against the truth
  0.990 and 0.991). It fails in `refine3D_states_cavgs`, whose states collapse
  into one; the unchanged code collapses the same way. Depending on the run,
  the failure shows as one state holding all sub-class averages, or as a state
  left with particles in one half only, which the distributed PCG solve
  refuses ("distributed PCG requires both halfsets"). Out of this run's
  acceptance by ruling.
- **The `solve3D_addon` pose-error excess** rose from 0.50° (Phase 2) to
  2.83° (Phase 3) and 4.76° / 3.78° (Phase 4), against a limit of 5°. Phase 4
  should not change it (by construction the streaming commit sums the same
  rows in the same order and writes the same spectra), and two
  runs of one build differ by 1°. Whether Phase 3 shifted it can be settled
  with the same repeated-run comparison as for `simulated_workflow_6vxx`.
- **Smaller items.**
  - The bootstrap master (`simple_commanders_euclid_distr.f90`, not in this
    phase's files) still writes its records through one open and sync per
    block; the new one-handle primitives would make it one sync.
  - The rebuild check now runs at every `refine2D` entry, including each
    `solve2D` stage and each stream pool round (a validation of the store's
    header and global section). A defect in the stream's store merges would
    show as repeated reseeding rather than a stop.
  - The generated indexes under `doc/code_overview/` still name the removed
    routines; they are regenerated by the maintainer.

## Merged into master

Master had moved on while the run worked. Two concurrent commits touched the
same code, and the merge into `~/src/SIMPLE` keeps them (the plan's section
15 has the detail):

- The sigma2 store carries no checksums any more, as master decided the same
  day: no header or global-spectra checksum, no deep validation, and the
  manifest identifies the committed sigma2 state by its generation. The
  streaming commit and the tests follow.
- `cls_expansion` keeps master's rewrite (conditional sub-class averages).
  From the run it takes the sigma2 interface (one consumability check, each
  member's weight from its half's spectrum), and the state-0 flags for
  sub-class halves without members, re-done in the rewritten code. The per-half
  weight table is not carried over.

One change was decided at the merge:

- The bootstrap's average image is now the mean of the images summed. The run
  preserved the old normalisation, which divided the sum of at most 25,000
  sampled images by the number of all active particles; above 25,000
  particles that shrank the subtracted power by the square of the sampled
  fraction (about 1/900 at 750,000 particles). Each worker now writes the
  number of images in its Fourier sum beside its per-half sums, and the master
  divides by their total. Up to 25,000 particles the seed changes only when a
  sampled spectrum was invalid and left out of the sum.
