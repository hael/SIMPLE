# Report of the run `speed2d3d_mac`: 2D tables, 3D table, NU filter and 2D look and feel

Report, 2026-10-19. The plan, with its progress per phase, is
`speed2d3d_run.md` beside this note. The run worked unattended on the
maintainer's Mac from 11:20 to about 18:30, on master `af13bc828`. Nothing
is committed; the changes are one diff (`full.diff` in the run directory
`~/agent_runs/runs/speed2d3d_mac`, per-item diffs in its `review/`).

## What was learned

- **The 2D tables can be smaller and read once, with results unchanged.**
  The dense 2D part file now stores 8 bytes per entry plus a short list of
  refined shifts (R2's part: 432 kB to 155 kB), and the master's dense table
  shrinks by the same factor. The sparse part file carries each particle's
  candidate count, so the master reads every candidate once. Both kinds of
  part are read in parallel, with a size check, a format word and a check
  that every particle is in exactly one part.
- **An uncached `solve2D` can equal a cached one, byte for byte.** Since
  `af13bc828` the 2D matcher prepared particles from the in-memory cropped
  image, but the 2D table still prepared them at the full box, so `cache=no`
  departed from `cache=yes` at the first table iteration. Now both use the
  same cropped image, and in shared memory the table hands its cropped
  particles to the matcher through the particle buffer: the matcher reads no
  particle in a table iteration. Uncached results change by rounding
  (accepted with the package).
- **The NU filter's label optimisation is about four times faster**, with
  the same label field. Reusing the setup's smoothings, skipping zero lines
  and planes, filtering the strided directions a plane at a time and
  walking packed voxel lists per colour took the label optimisation of the
  three-state case from 23.6 s to 5.5 s (one thread, an indication).
- **The shift matrix cannot be fused into the alignment kernels identically
  on this compiler.** The fused kernels compute every factor with the same
  statements as the shift matrix, yet results differ (one SW1 particle's
  shift in the 7th digit; E2, E2S and R2 broadly): the compiler contracts
  multiply-adds differently in the two contexts. The item is restored.
- **The 2D log now reads like the 3D log**: one `>>> CONVERGED: yes|no` per
  iteration with the decision that actually stops the run, the matching
  low-pass, stage lines that are not cut, entering and returning lines, one
  line per pre-step. E2's log shrank from 872 to 760 lines.
- **The test harness can lie quietly.** The first item build was synced
  with file times older than its objects, so make did not recompile it, and
  every case was "identical". The new section of the item's bench report
  was missing, which gave it away. Items are now synced by content, and the
  run checked each item's observable trace in its outputs.

## How the run tested

- **Identity**: every case at one thread with a fixed seed, in fixed
  directories, all regular files compared except the exclusions declared in
  Phase 0 (`exclusions.tsv`: project files through their dumps, benchmark
  reports, process ids and times, the solve3D manifest's run id, PCG
  timing lines). Table part files are compared through a lossless decoder
  that knows every layout, so a layout change is judged by the
  information it carries.
- **Cases**, all made on the base build and each run twice there without a
  difference: the simulated workflows SW1 (614 files) and SW6 (546); the
  three-state `refine3D` MSI (216), the correlation objective with
  fractional update CCI (50), the PCG backend PCI (53), the NU programs NU3D
  (14); `solve2D` on 1,500 simulated particles, two parts E2 (140) and
  shared memory E2S (137); `refine2D refine=prob` R2 (30); E2 and E2S with
  `cache=yes` (E2C, E2SC); AD1 (46 check values) and CR1 (60).
- **Two case decisions.** MSI, CCI and PCI take `box_crop=128`: without it
  `refine3D` does not crop and the cropped-record path of `af13bc828` is
  never run. E2CN and E2SCN run E2 and E2S with `cache=no` in the
  directories of E2C and E2SC, because `metadata.json` and the project dumps
  hold the case path.
- **`af13bc828` checked.** Its parent `aea5b1846` and the base gave
  byte-identical MSI, CCI and PCI (the 3D record path, uncompiled until the
  run).
- **Parallel code**: every new or changed parallel region has a tester check
  at the Mac's 11 threads, repeated 20 times or more, bit for bit where the
  result is determined (the 2D part reads, the sparse 2D assignment, the
  tent smoothing; the 3D read already had one), and the label field at a
  fixed thread count. The cropped 2D batch ran in E2S at 11 threads 20
  times, each to the end with finite class averages.
- **Phase ends**: SW6, `./compile_debug.sh` on an idle Mac (the fast gate
  passed every time, 14.6 to 16.0 s), and E2S at 11 threads and MSI with two
  parts of 5 threads, which always ran to the end with finite results.

## What changed

### The 2D tables and the master's work between them (Phase 1; kept, identical)

- **Timers.** The `PROB_ALIGN2D` report gains a section splitting the reads
  (the parts' heads, their candidates) and the assignment (class lists:
  counts, fill, sort; the frontier loop; the leftover particles). On E2 every
  phase is 0 to 3 ms of a 2 s pre-step.
- **One read of each sparse part.** The 2D part file (format word -4)
  carries the candidate count of every particle; the master sizes its
  store from the counts and reads every candidate once.
- **The compact dense record.** The dense 2D part file (format -5) holds one
  row of (distance, in-plane index) per particle and class and the refined
  entries (particle, class, shift). The master keeps 8-byte rows and the
  refined entries in compressed-row form, and rebuilds a full candidate when
  it assigns. The 2D tester writes its part files through the table's own
  writer and checks that a refined shift survives.
- **Parallel part reads**, both layouts: an OpenMP loop over the parts with a
  status per part raised after it, file names built before it, a size check
  against the header, a format word, and an exactly-one-part check with
  `atomic capture`. Tester: seven parts read at one thread and 20 times at
  11 threads, equal bit for bit.
- **The sparse assignment's class lists** are gathered with their particles
  and distances during the fill, sorted in parallel over the classes with
  the same unstable sort, and walked in order; peak counts are computed
  once. The frontier loop is unchanged. The tester keeps the replaced
  algorithm and matches it on tied tables at four seeds, and at one against
  11 threads.
- **Not done**: the top-k sampler of the frontier loop. Its condition (the
  loop dominating) was not met on E2; the new timers settle it on a large
  case.

### 2D look and feel (Phase 2; kept, outputs unchanged)

- **Log lines.** `check_conv2D` tests the criteria; the 2D matcher (and the
  stream's 2D pool through it) and the distributed strategy each print one
  verdict after their own `minits`, `maxits` and coverage rules. The
  meaningless `PROB NEIGH MODE` line and the spread of a constant low-pass
  are gone; the dead `max_cls_pop` and its stage-table column are deleted.
  The stats file's keys are unchanged.
- **The solve2D interface.** Correct help for `center` (`yes`), `refine` (how
  the three modes differ), `lpstart`, `lpstop` and `cenlp` (the mask-diameter
  rules); `nsample` is Standard and shows 200,000; `nstages` is exposed with
  its 2D meaning; the unused `refine2D` program variable is deleted.
  `default_audit.py` follows `exec_solve2D` into `exec_solve2D_workflow`, so
  the interface shows the defaults the code applies (`autoscale=yes`,
  `cls_init=rand`). Only solve2D's entry of the UI JSON changed (compared
  before and after).
- **Class-average names** come from `src/defs/simple_refine2D_fnames.f90`
  at every non-stream site, including the parse of the iteration in
  `make_cavgs`; the dead `refine3D_iter_refs_fname` is deleted.
- **Shared sigma2 helpers.** `simple_strategy_sigma2.f90` holds the store
  preparation and the bench counts both strategies copied.
  `ensure_canonical_sigma2` stays per strategy: its shared core is already
  one call, and the rest is what differs.
- **`solve2D_policy.md`** has Defaults and Stage-controller sections, the
  canonical sigma2 store, and the fill-in as the code does it.

### The 3D table and the NU filter (Phase 3)

- **The 3D master table by reference (kept, identical).** `dists` and
  `inpls` are (particles, references), so each reference's distances are
  contiguous; a part's rows are read in chunks of up to 8 MiB and scattered;
  refined entries are counted during the parallel read and placed in
  parallel; the state labels and the invalid/flat check run references
  outside with per-particle best-so-far arrays. The dead `score_min`,
  `score_spread` and `score_mode` are gone.
- **The shift matrix fused into the alignment kernels (not kept).** Restored
  as the plan says; the dead store in the correlation gradient went with it.
- **The NU label optimisation (kept, identical).** Sub-timers (level
  smoothings, Potts sweeps) in `VOLASSEMBLE`; the setup's smoothed candidate
  and its smoothing norm reused at each level (42 smoothings to 28); the
  tent smoothing skips zero rows, planes and slabs and filters Y and Z a
  plane at a time (tester: bit for bit against the old filter, at one thread
  and 20 times at 11); packed voxel lists per colour in the Potts sweeps
  (tester: the label field 20 times at 11 threads). The four parts were
  judged together, because their own snapshots held the kernel change.

### 2D cropped particles shared (Phase 4; kept, judged against `cache=yes`)

- `crop_particles2D` decides for the 2D table and the 2D matcher. The table
  prepares alignment from the in-memory cache entry, as a cached particle.
  In shared memory it opens the particle buffer for records instead of raw
  images; the matcher copies a batch's records and neither reads nor
  prepares the full box. E2CN equals E2C (140 files) and E2SCN equals E2SC
  (137 files) byte for byte. In distributed runs a table worker crops each
  particle itself: one more inverse and forward transform at `box_crop` per
  particle than the cached path.

## Time, as an indication

One thread unless said; the Mac was shared with other jobs, so differences
of a few seconds are noise.

| case | base | final | note |
| --- | --- | --- | --- |
| MSI, step 1 + step 2 | 59 + 55 s | 43 + 42 s | the NU label optimisation over all volume assemblies 23.6 s to 5.5 s |
| PCI | 42 s | 40 s | |
| E2, `cache=no` (E2CN) | 84 s | 80 s | |
| E2S, `cache=no` (E2SCN) | 121 s | 111 s | the matcher's particle preparation 20.8 s to 15.9 s; reads 0.8 s to 0.3 s |
| E2S, 11 threads | 24 s | 23 s | |
| MSI, two parts of 5 threads | 22 + 23 s | 18 + 23 s | |

The 2D table and assignment timers matter on large cases (msp1: 141,838
particles, 150 classes), not on E2.

## Not kept, not done, and questions for the maintainer

- **Not kept**: the shift matrix fused into the four alignment kernels (not
  bit-identical, above).
- **Not done**: the top-k frontier sampler (its condition was not met).
- **Question: derived defaults in the interface.** The numeric defaults shown
  for solve2D's `lpstop` (8) and `cenlp` (20) are parameter baselines that
  the 2D code never applies; it derives them from the mask diameter. An
  empty numeric default is serialised as 0.0, and a GUI that sends shown
  defaults would then run `lpstop=0`, so the shown values were left and only
  the help gives the rule. The UI schema needs a "derived" default.
- **Stream code (not edited)** builds five class-average names inline that
  `simple_refine2D_fnames.f90` could build: `simple_stream_pool2D.f90` (the
  JPEG name, twice; the references, twice) and
  `simple_stream_refine2D_utils.f90` (the prefix of a pool iteration). The
  stream pool prints the new 2D convergence verdict through the matcher.
- **The NU unit suite** takes about one second longer on an idle Mac (the
  label-field repeat test).

## Proposed skill text

- `simple-solve2d`: "The 2D part files have a format word: -4 for the
  candidate stream (with per-particle counts), -5 for the dense rows (8 bytes
  per entry plus the refined entries). `prob_align2D` reads them in parallel
  and checks size, format and that each particle is in exactly one part.
  When `crop_particles2D` holds (uncached, not the stream pool, `box_crop <
  box`, ptcl2D), the 2D table and the matcher prepare particles from the
  in-memory cache entry, and an uncached run equals a cached one byte for
  byte; in shared memory the table's records feed the matcher."
- `simple-solve2d`: "`check_conv2D` does not print a verdict: the caller
  prints `convergence%print_verdict` after its own minits, maxits and
  coverage decision."
- `simple-refine3d`: "The master's dense table is `dists`/`inpls(nptcls,
  nrefs)`, each reference's distances contiguous."
- `simple-main-nu-filt`: "The level loop of `optimize_nu_cutoff_finds` reuses
  `dmats_mask` and `norms_mask` for the level's own label; the tent smoothing
  is checked bit for bit against a column-by-column reference in the
  tester."

## Files

Changed: `simple_eul_prob_tab2D.f90` and its tester, `simple_eul_prob_tab.f90`
and its tester, `simple_commanders_prob.f90`, `simple_convergence.f90`,
`simple_strategy2D_matcher.f90`, `simple_matcher_ptcl_batch.f90`,
`simple_matcher_ptcl_io.f90`, `simple_refine2D_strategy.f90`,
`simple_refine3D_strategy.f90`, `simple_commanders_solve2D.f90`,
`simple_solve2D_controller.f90`, `simple_ui_refine2D.f90`,
`scripts/default_audit.py`, `simple_refine3D_fnames.f90`,
`simple_commanders_mkcavgs.f90`, `simple_starproject.f90`,
`single_commanders_nano2D.f90`, `simple_nu_filter{,_bank,_state,_potts,_vars,_tester}.f90`,
`simple_nu_state_filter.f90`, `simple_commanders_rec_distr.f90`,
`simple_tent_smooth.f90`; new: `src/defs/simple_refine2D_fnames.f90`,
`src/main/strategies/parallelization/simple_strategy_sigma2.f90`. Documents:
`solve2D_policy.md`, `particle_cache_policy.md`, `refine3D_policy.md`,
`nonuniform_filtering_policy.md`, `particle_read_path_policy.md`,
`speed_architecture_proposals.md`, the plan and this report.

Evidence in `~/speed2d3d_mac` (logs/iNN_summary.txt, p1-p4_end.txt,
i21_judge.txt, i21_stress.txt; outputs in out/base and out/item) and in the
run directory (journal.md, exclusions.tsv, phase0_determinism.txt,
review/).

## Commit message

```
2D tables, 3D table, NU label optimisation, 2D log and interface

2D probability tables: per-particle counts and a format word in the part
files, a compact dense record (8 bytes per entry plus refined shifts),
parallel part reads with size, format and coverage checks, the sparse
assignment's class lists sorted in parallel, read and assignment timers.
The 2D table prepares particles from the matcher's in-memory cropped image,
and in shared memory hands them to the matcher through the particle buffer,
so cache=no equals cache=yes. 3D master table stored by reference. NU label
optimisation: setup smoothings reused, tent smoothing skips zeros and
filters planes, packed Potts sweeps. 2D log lines and convergence verdict
as in 3D, solve2D interface texts, class-average names through helpers,
shared sigma2 helpers, policy updates. Identical at one thread except the
uncached 2D path, which now equals the cached one.
```
