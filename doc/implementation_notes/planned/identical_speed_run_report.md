# Identical speed changes: report of the run `speed_mac`

Report, 2026-10-17; Phase 6 added 2026-10-18. The run carried out the plan in `identical_speed_run.md`
on the maintainer's Mac, unattended, as one uncommitted diff. "Identical"
means: the outputs of the identity cases are byte for byte those of the
build before the run.

## What we learned

- **Every planned item is identical, and all but one are kept whole.**
  Fourteen changes in four areas (the probability tables, the polar
  kernels, the nonuniform filter with volume assembly, the master's project
  handling) leave every output of the one-thread identity cases unchanged.
  The fourteenth is kept in part (below).
- **The gains are small on these small cases.** The identity cases have 87
  particles and box 192. Under the same conditions on the Mac, their wall
  times are within about 2% before and after; the neighbourhood table's
  refinement is about 9% faster, the filter about 1%. The changes target
  large runs: big part files, many references per particle, many particles
  in the project. The gains must be measured on the test machine (last
  section).
- **Two deciding timings.** Filtering the candidates again in the filter's
  composition is faster than writing and reading them (9.9 s against
  10.5 s, box 300). Copying a 751,729-particle project in memory takes
  0.28 s, reading it 1.3 s.
- **The identity cases caught three of my own mistakes before any item was
  judged.** A timer array freed before its report was written (a crash in
  the three-state case), killed volumes reused by the kept reprojection
  builder (a crash in the second iteration), and that builder writing the
  particle field into the wrong project file: the shared input project of
  four cases. The project was restored byte for byte from an identical run
  (its checksum matches the original), the runs that had read it were
  discarded, and every case now checks that input after it runs.
- **One-thread determinism held throughout.** Every case gave the same
  bytes twice on the base build, and every later comparison was exact
  under the exclusions declared at the start. No case was rerun until it
  passed; no tolerance was used.

## What changed

### The probability tables

- **The neighbourhood table read in one pass.** A table part now writes
  the number of candidates of each particle after the particle indices.
  The master sizes its store from these counts instead of scanning every
  candidate, then places the candidates in file order as before (the order
  decides ties).
- **A compact record for the dense table** (the `prob` table, where every
  particle has a value for every reference). Each particle has one row of
  (distance, in-plane index) per reference, 8 bytes, and a short list of
  the references whose shift was refined. The file says which layout it
  has (a format word in the header). The worker writes a batch's rows at
  the particles' positions and checks at the end that every row is there.
  The master holds the rows and the refined lists and rebuilds the entry
  the assignment reads (`dense_candidate`). The master's table is no
  longer filled twice before the read. On one CCI part file: 532 KB before,
  153 KB after.
- **The dense part files read in parallel.** Two new checks follow the
  loop: every particle of the master is in exactly one part, and the parts
  agree on the rotation grid. Before, a missing particle silently kept
  default values.
- The old candidate writer serves only the neighbourhood table now and
  moved into its type; each empty part is written by its own table type.
  The unused `fill_tab` of both types is gone.

### The polar kernels (the shift search and the score loop)

- **The shift matrix converted only on the searched shells.** Every
  consumer reads only those shells.
- **The repeated objective reused.** After L-BFGS-B, the shift search
  evaluated the objective over all rotations at the point the minimiser
  ended on; the minimiser's callback had just evaluated the same vector at
  that point. Inside `minimize` the last vector is kept with its shift and
  reused when the shift is bitwise equal.
- **One double-precision copy of the reference per minimisation**, in a
  per-thread buffer, instead of one per evaluation. The unreachable
  gradient-only path and `shift_ref`, which had no caller, are deleted.
- **The score loop without allocations.** The compiler confirmed both
  suspicions of the plan: a heap trampoline per call for the internal
  distance functions and an allocation of the score vector per call. The
  likelihood routes now use per-thread buffers and a sampler over an array
  of distances, which shares its draw with the old sampler (still used by
  the assignment).

### Volume assembly and the nonuniform filter

- **Timers.** The volume-assembly report has a new section: the filter
  split into candidate filtering, objective, label optimisation and
  composition, and per state the restore and the filter. The summary
  script shows them. On MSI's three-state step the label optimisation is
  the largest part (8.7 of 14.9 s).
- **Fewer reads.** Two reads of the restored halves that nothing used are
  gone, and the filter takes its input halves from memory instead of
  reading back the files just written.
- **The objective right after each candidate** is filtered, without
  writing and reading the candidate.
- **The composition without scratch files.** The filtered references are
  made by filtering the input pair again per label, with the same steps as
  the candidates, so nothing writes the filter's scratch files any more.
  The unused `nu_filter_vol` is deleted.

### The master's project handling (distributed `refine3D`)

- **Fewer rereads.** Volume assembly no longer rereads the particle field
  for its trailing counts (Phase 6 removed the last case, below).
  The matching low-pass is read from one particle record instead of the
  whole segment; a record range given to `read_segment` is ignored for
  project files, so this goes through the binary reader directly.
- **The reprojection builder kept for the run.** Its toolboxes are made
  once; later iterations read only the project segments the model uses and
  release them afterwards.
- **`prob_align` gets a copy of the master's project.** It runs right after
  the master's full read, so the copy is what it would read.
- **Volume assembly takes the master's project and hands its own on**
  (Phase 6, below).

## How identity was shown

- **The cases** (all at one thread, fixed seed, the CTest environment):
  the simulated workflow of beta-galactosidase (SW1, both its
  shared-memory and its two-part step), a three-state `refine3D` with the
  dense table and then the neighbourhood table (MSI), correlation-objective
  refinement with fractional update and trailing reconstruction (CCI), the
  PCG backend with the filter (PCI), and the filter's standalone programs
  (NU3D). At each phase's end also SW6 (the spike), `solve3D_addon` (AD1)
  and `cont_refine3D_1jxy` (CR1).
- **The comparison:** every regular file of a case's output, byte for
  byte. Declared once, at the start, with reasons
  (`exclusions.tsv` in the run directory): project files (compared through
  text dumps of their segments), benchmark reports, the time and process-id
  fields of the GUI metadata, job records holding process ids, previews
  named by process id, the run identifier of the `solve3D` manifest, and
  the timing lines of the PCG diagnostics. AD1 and CR1 work in directories
  named by process id; their check values are compared (44 and 59).
- **Each item was built from its own snapshot** and its five cases run
  against the base outputs. Phase-end cases ran on the phase's last build:
  SW6, AD1 and CR1 were identical after each of the four phases. Every
  comparison came out identical.
- **A decision for the maintainer to review:** the table part files stay in
  the case directories, and the table items change their layout on
  purpose. They were compared through a lossless decoder (every field, bit
  for bit, particle by particle; only the points where thread buffers were
  flushed are not kept) instead of byte for byte. No exclusion was added.
- **Exact checks in the testers** (the fast gate's unit suites): the dense
  part files read back exactly; the parallel read at one and at four
  threads; the scores independent of the shells above the searched band;
  the loaded reference evaluating exactly as the copied one; the shift
  searches and likelihood scores at one and at four threads; the array
  sampler drawing exactly as the old one; every composed voxel of the
  filter equal to its label's candidate. The fast gate passed after every
  item and at every phase's end (15 of 15 suites, 12.7-14.1 s).

## Timings

The two timings that decided an item, each with nothing else running:

| what | before | after |
| --- | --- | --- |
| `nu_filt3D`, PfCRT half maps, box 300, 8 threads, three runs interleaved | 10.56, 10.58, 10.27 s (scratch files) | 9.95, 9.85, 9.75 s (filtering again) |
| the largest project at hand, 751,729 particles, 438 MB, three times | read 1.43, 1.28, 1.27 s | copy 0.31, 0.28, 0.25 s |

Peak memory was the same before and after the filter change (1.89 GB).

One-thread benchmark lines, as an indication only (MSI, the Phase 1 build
against the final build, the same load on the Mac): the dense table's fill
8.64 against 8.49 s, the neighbourhood table's fill 19.30 against 17.79 s
(its refinement 14.43 against 13.20 s), volume assembly's filter 102.6
against 101.3 s. A comparison with the base runs shows larger gains, but
those runs shared the Mac with seven to nine other jobs and are not a fair
reference.

## Not done, and why

- **The states' volume assembly in parallel** (from the week plan) is not
  an identical change and was not attempted: several sums of a state's
  assembly depend on the thread count, the filter is module state, and the
  policy on multi-state peak memory forbids holding several states'
  volumes. Per-state timers were added instead.
- **Two heap trampolines remain**, in the bounded and power samplers used
  by the 2D sampling strategies and `greedy_sub`, outside the table fill.
  The array sampler of this run would remove them.

## Decisions for the maintainer

All four were ruled on 2026-10-18 (Phase 6, below).

1. Accept the decoded comparison of the table part files (above).
   Accepted.
2. Volume assembly rereads the particle field when it has back-filled
   missing projection labels, and that reread throws the labels away; the
   field is written later. So a trailing run loses the labels and others
   keep them. Fixed in Phase 6.
3. Whether volume assembly should hand its project to `commit_sigma2` and
   the master (the remaining project reads). Done in Phase 6.
4. On this Mac, AD1 fails its own floor at one thread
   (`cohort_pose_error_excess_deg` 6.00 against 5.0) on the base build and
   on every build after, identically; the run `ms_speed` found that floor
   sensitive to the draw. SW6 fails its volume-correlation floor (0.780
   against 0.80), as the plan expected. SW1 and CR1 pass. These floors
   are left to the run `test_truth` on the Dell.

## Phase 6: the rulings of 2026-10-18

The maintainer accepted the decoded comparison of the part files, left the
Mac's failing floors (AD1 at one thread, SW6) to the run `test_truth`, and
asked for two more items, as fast as possible.

**Only the unit tests were run.** By ruling, the gate was the Debug build
and the fast-gate unit suites (15 of 15 pass after each item, and after
the final `./compile_debug.sh`). No identity case, Release build or timing
was run. Each item has its own snapshot and diff
(`review/item_21_assembly_hands_on_project.diff`,
`review/item_22_assembly_keeps_backfilled_labels.diff`), so the Dell can
test them apart.

- **Item 21: volume assembly hands its project on (ruling 3).** In the
  distributed gridding path (not `refine=sigma`, which skips the merge),
  the master passes its project to volume assembly (`execute_hand_on`),
  which copies it instead of reading the file. The master read it at the
  start of the iteration, `merge_algndocs` wrote the field from it, and
  nothing else writes the project before assembly, so it equals the file.
  Assembly hands its own project back as it ends. The low-pass refresh
  reads the one record it needs from it, and `commit_sigma2` takes a copy
  (`execute_from_project`) instead of its full read. After assembly the
  master writes only the `out` segment, which `commit_sigma2` does not
  use. No write of the project changed. The shared-memory path is
  unchanged. Meant to be identical.
- **Item 22: the back-filled labels are kept (ruling 2).** The trailing
  counts use the field in memory; the reread is gone. The back-filled
  projection labels do not enter the counts, so the counts are the same,
  and the field volume assembly writes now keeps the labels on every path.
  A fix: results change where labels are back-filled.

**What the Dell must check.**

- Item 21: byte identity against the build before it on SW1, MSI, CCI and
  PCI at one thread.
- Item 22: which callers reach the back-fill (volume assembly on a project
  whose particle field has no `proj` labels, on a trailing run:
  `trail_rec=yes`, or a `trail_seed` run); that the only differences are
  the kept `proj` labels in the written field and what follows from them.

## Files

Changed: `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`,
`simple_eul_prob_tab_neigh.f90`, `simple_eul_prob_tab_utils.f90`,
`simple_eul_prob_tab_tester.f90`;
`src/main/commanders/simple/simple_commanders_prob.f90`,
`simple_commanders_rec_distr.f90`, `simple_commanders_resolest.f90`;
`src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`,
`simple_polarft_corr.f90`, `simple_polarft_geom.f90`,
`simple_pftc_shsrch_grad.f90`, `simple_pftc_inplane_tester.f90`;
`src/main/nu_filt/simple_nu_filter.f90`, `_vars.f90`, `_state.f90`,
`_bank.f90`, `_apply.f90`, `_tester.f90`;
`src/main/volume/simple_nu_state_filter.f90`; `src/main/simple_builder.f90`;
`src/main/strategies/parallelization/simple_refine3D_strategy.f90`;
`scripts/plot_refine3d_bench.py`; the two policies and three notes in
`doc/`. None of the files the plan reserves for the run `test_truth` was
touched. Files that run may also edit: `simple_eul_prob_tab.f90`
(`ref_assign` reads the distances from the new rows; its algorithm is
unchanged), `simple_commanders_rec_distr.f90` and
`simple_refine3D_strategy.f90`.

Proposed skill text (not applied): in `simple-main-nu-filt`, replace
"Always preserve cleanup of module allocatables and scratch cache files"
by "Always preserve cleanup of module allocatables; the filter writes no
scratch files (the composition filters the source pair again)". In
`simple-refine3d`, under the probability tables: "The dense table's part
file is a row of (distance, in-plane index) per particle and reference plus
a refined-shift list; the neighbourhood and 2D tables keep candidate
streams."

## Timing cases for the test machine

First identity again on the x86 build, which this Mac's ARM build does not
prove: SW1, MSI and CCI at one thread on the base and the final build.
Then, interleaved, twice each way, the cases of the runs `prob_speed` and
`ms_speed`: A and A-shm (`solve3D` of beta-galactosidase), MS (three
states from the class averages), M (the dense three-state table at 1,000
directions) and E (`solve2D`). Read in their benchmark reports: the table
fills (score loop and refinement), the pre-step's table reads, the
volume-assembly filter and its new split, and the strategy's
reprojection-model time.

## Evidence

The run directory `~/agent_runs/runs/speed_mac`: `journal.md` (every step,
decision and defect), `review/item_NN_*.diff` (one diff per item),
`full.diff`, `exclusions.tsv`, `phase0_resources.txt`,
`phase0_determinism.txt`. Scratch `~/speed_mac`: `out/<build>/<case>_<label>/compare.txt`
for every comparison, the timing data in `timing/`, the helper scripts in
`bin/`.

## Commit message

    Identical speed changes to the tables, kernels, filter and master

    Compact dense part files read in parallel, a one-pass neighbourhood
    read, one reference copy per shift search with the repeated objective
    reused, no allocations in the score loop, the nonuniform filter without
    scratch files and with timers, fewer project reads by the distributed
    master: volume assembly takes the master's project and hands its own
    to the low-pass refresh and commit_sigma2. Byte-identical on the
    one-thread identity cases, except the last hand-off, so far checked by
    the unit suites only.

    Fix: volume assembly keeps the projection labels it back-fills; the
    trailing counts no longer reread the particle field, which dropped
    them on trailing runs.
