# 2D tables, 3D tables and the NU filter, and 2D look and feel: the Mac run `speed2d3d_mac`

Implementation note, 2026-10-19. **Carried out 2026-10-19** (progress below;
report in `speed2d3d_run_report.md`). The execution
plan of the unattended ten-hour run `speed2d3d_mac` on the maintainer's Mac.
The design is in `speed_architecture_proposals.md` (sections 1.2 b and d, and
1.3) under `doc/implementation_notes/planned/`; what the previous Mac run did
is in `doc/implementation_notes/completed/speed2D_run_report.md`. This note
holds the run's items, tests, phases and progress. Where this note and the
proposals note disagree, this note holds: it was checked against the code at
`af13bc828` on 2026-10-19 (section "What the reading found").

## Rulings of 2026-10-19

- The chunk: finish the 2D work (the 2D tables and the master's work between
  them, then the 2D look and feel), and the 3D table, shift-search kernel and
  NU filter items. Every item of Phases 1 to 3 leaves results byte-identical
  at one thread.
- Last, only if time is left (Phase 4): the 2D table prepares its particles
  from the same in-memory cropped image the 2D search has used since
  `af13bc828`, and the shared-memory particle buffer keeps those cropped
  images. This changes uncached results by rounding (accepted with the
  package); its identity is judged against `cache=yes` runs of the base
  build, which must come out byte-identical.
- Testing as in the previous run: one-thread identity on the simulated cases,
  and for every new or changed parallel region a thread stress test at the
  Mac's full thread count, repeated; one multi-threaded run per phase to catch
  crashes and races.
- Not in this run: the master drawing the stage 1-2 sample, partition
  balancing by active particles, the `maxits`/`minits`, `fillin` and
  `prob_snhc` renames, the particle-cache and `nsample` defaults, "one project
  owner on the master", the correlation's shift-independent term (it is not
  identical, below), and any edit of stream code (`src/main/stream`,
  `commanders/stream`: Joe's area; list what would change there in the
  report).

## What the reading found

Checked on 2026-10-19 at `af13bc828`; the run trusts these only after
re-reading the code.

- **The master's 2D pre-step time is the sparse table.** On msp1 (141,838
  particles, 150 classes, 10 parts) the 3.8-4.3 s of `ref_assign` and the
  1.6-1.8 s read per iteration fall in the `prob_snhc` iterations (solve2D
  stages 3 to 5), which use the sparse 2D table; the dense stage 6 assigns in
  1.09 s. The parts score 131-148 of 150 classes per particle, so the sparse
  table is nearly full (`~/Processing/msp1/BENCH/PROB_ALIGN2D_BENCH.txt`,
  read-only).
- **The 2D tables.** The master holds `loc_tab(nclasses,nptcls)` of
  `prob_candidate` (six 4-byte fields, 24 bytes; 511 MB on msp1), filled in
  parallel inside the timed read. A part file is a header, the particle
  indices, the seed table, then chunks of particle indices and 24-byte
  candidates (28 bytes per entry). The dense reader (`read_tab_to_glob`) runs
  part by part, serially, with no file-size check, no format word and no
  check that a particle is in exactly one part. The sparse reader
  (`read_tabs_to_glob` then `read_sparse_tab_to_glob`) reads every part twice:
  once to count candidates per particle, once to read them. Reading arrays of
  a derived type may move one component at a time; time it before relying on
  it.
- **The sparse 2D assignment** (`ref_assign_sparse_likelihood`): building the
  class lists is serial (a count pass and a fill pass over all candidates,
  then per class a gather of distances and an `hpsort` of about 140,000
  values); the frontier loop then runs once per particle (`sample_likelihood_index`
  sorts about 150 values and draws one `ran3()`), scans the particle's
  candidates three times (`search_frac`, `search_npeaks`,
  `update_frontier_after_assignment`) and walks class heads with random reads.
  The only random draw is one `ran3()` per loop iteration, in assignment
  order. `hpsort` is not stable: its order among equal distances decides a
  class's head, so it stays. The swap-remove in `sync_frontier_class` decides
  ties within a particle, so the order of classes within a particle stays.
- **The 3D master table is already particle-major.** `rows(nrefs,nptcls)`
  keeps each particle's references contiguous. The assignment's main loop
  sorts each reference over all particles and reads `rows(iref,i)%dist` with
  a stride of `nrefs*8` bytes; that is the loop that does not stream. The
  useful change is the opposite layout (each reference's distances
  contiguous). `read_dense_part` issues one read per particle, though a
  part's rows are one contiguous block in local order.
- **The correlation's shift-independent term is not identical when hoisted.**
  The reference norm in the correlation kernels is computed from
  `pft_ref*shmat_8*ctf`; the shift matrix's modulus drifts by rounding along
  its recurrence and the term depends on `irot`, which changes during an
  alternating search. Computing it once per minimisation changes results.
  Only fusing the shift matrix into the kernel loops can be identical.
- **The NU filter's label optimisation** (`optimize_nu_cutoff_finds`) runs 35
  tent smoothings at 8 candidates plus up to 7 re-smoothings of the mask norm;
  the level-1 values equal the setup's smoothed candidates bit for bit, so
  reuse cuts 42 smoothings to 28. Skipping all-zero lines and blocking the
  strided passes of the tent smoothing, and packed voxel lists per colour in
  the Potts sweeps, keep every operation; replacing the neighbour sum by
  label-histogram counts would not. `estimate_nu_label_smooth_beta` uses an
  OpenMP real reduction, so the filter is identical only at a fixed thread
  count (as today).
- **2D look and feel.** `check_conv2D` prints `CONVERGED: .YES.`/`.NO.` before
  its callers apply `minits`, `maxits` and full-assignment coverage, so the
  printed verdict can contradict whether the run stops; it prints a `PROB
  NEIGH MODE` line that means nothing in 2D and the SDEV/MIN/MAX of a
  low-pass that is the same for every particle. The solve2D stage line's
  format `(A,I3,A9,F5.1)` cuts `' WITH LP ='` to `>>> STAGE   3 WITH LP  12.5`.
  No line marks entering and returning from a `refine2D` stage. `prob_align2D`
  and the 2D assignment print 7 to 10 lines per iteration. The stage table's
  `MAXCLSPOP` column prints the dead `max_cls_pop` (always 0). On the
  distributed master `params%lp` is never updated (the workers set the band),
  so a 3D-style "this iteration" low-pass line must print the particles'
  `lp` field. `check_conv2D` also writes the stats file's keys (an output:
  unchanged) and is called by the stream pool (its log lines change there
  too; nothing parses them).
- **The solve2D interface** shows defaults that the code does not apply:
  `center` shows `no` (the code uses `yes`), `lpstart` 30, `lpstop` 6 and
  `cenlp` 30 (the code derives them from the mask diameter). The displayed
  numeric defaults come from `scripts/default_audit.py` at build time, which
  does not follow `exec_solve2D` into `exec_solve2D_workflow`, so no solve2D
  block is generated. `nsample` is Advanced in solve2D (Standard in solve3D)
  and shows 0 (the 2D default is 200,000). `nstages` is read but not exposed;
  in 2D it re-spaces the low-pass ladder (in 3D it truncates the plan), so the
  wording differs. `refine=snhc_smpl` is listed without explanation: it
  differs from the default only in the last stage. `type(ui_program) ::
  refine2D` in `simple_ui_refine2D.f90` is never constructed (dead).
- **Class-average file names** are built inline in about fifteen places (the
  solve2D commander and controller, the refine2D strategy, the 2D matcher,
  `make_cavgs`, the star-file writer, nano2D; stream has five more), some
  with `params%ext` and some with `MRC_EXT`; `make_cavgs` reads the iteration
  back out of a name by offset. `simple_refine3D_fnames.f90` is the
  precedent; its `refine3D_iter_refs_fname` is called nowhere.
- **Sigma2 helpers.** `prepare_canonical_sigma_update` is identical in the 2D
  and 3D strategies; `record_sigma2_bench_context` differs only in the bench
  type; `ensure_canonical_sigma2` differs (objective check, iteration, label,
  a reuse line). Further variants live in the `make_cavgs`, `rec3D` and
  solve2D code.

## Why thread stress testing

The run `speed_mac`'s parallel read of the dense table parts passed one-thread
identity, yet held a race the maintainer fixed in commit 43de9b340: string
temporaries built inside the OpenMP loop and handed to `fopen`, a coverage
check that needed `atomic capture`, and a work array sized by `nthr_glob`
while its loop ran on `omp_get_max_threads()` threads. One thread runs an
OpenMP loop serially, so identity can never show thread safety. Rules for
every parallel region this run adds or changes:

- No string or allocatable-character temporaries, and no calls into routines
  with module-level buffers, inside the region: build file names and messages
  before it.
- Work arrays sized by the team the region actually uses.
- Errors inside the region set a status per task, raised after the region.
- A tester check that runs the routine at one thread and at the Mac's full
  thread count (`omp_get_num_procs()`), at least 20 times, and compares bit
  for bit where the result is deterministic by construction; where it is not
  (dynamic schedules that move random draws), every run must complete with
  finite results.

## The clock

The run has ten hours from its launch: `launch.sh` writes the launch time and
the deadline into the run's `state/` directory and each phase prompt states
them. The report and `full.diff` (Phase 5) are written before the deadline.

- Check the clock before starting each item. Do not start an item that
  cannot be finished, tested and snapshotted, and still leave the phase-end
  runs (about 25 minutes) and the 45-minute close before the deadline; write
  it up as not done, for time. In Phase 3 leave another hour for AD1 and CR1.
- AD1 and CR1 run at the end of Phase 3. If Phase 3 is cut or skipped, they
  run in the background at the start of Phase 5 and the report states their
  verdict (or that they did not finish).
- Phase 4 starts only if at least two and a half hours remain. The phase
  loop also skips a phase that would start less than 45 minutes before the
  deadline, so that Phase 5 always runs.
- Inside a phase the items are in priority order: when time runs short, the
  later ones go.

## The items

### Phase 1: the 2D tables and the master's work between them

1. **Timers first** (benchmark output only). Split the `PROB_ALIGN2D` report's
   read and assignment lines: the sparse read's two passes; the assignment's
   class-list building (counting, fill, per-class gather and sort), frontier
   loop and leftover particles, for both the sparse and the dense
   assignment. Run E2 and R2 on the timer build and record where the time
   goes on these small cases, as an indication only.
2. **One pass over each sparse part file.** The 2D part file carries the
   candidate count of every particle after the particle indices, as the 3D
   neighbourhood table does; the master sizes its ragged store from the
   counts and reads every part once. Read into intrinsic arrays (indices,
   distances, in-plane indices, shifts) and then build the store, if the
   timer shows the derived-type transfer is the cost.
3. **Parallel 2D part reads** for the sparse and the dense table, as 3D's
   `read_dense_part` does: an OpenMP loop over the parts with a status per
   part raised after the loop, a file-size check against the header, a
   format word, and a check that every particle is in exactly one part (with
   `atomic capture`). File names built before the loop. Stress test.
4. **The sparse 2D assignment's class lists in parallel.** The per-class
   gather and `hpsort` run in parallel over classes with work arrays sized by
   the team; gather the distances class by class during the fill pass; keep
   each class's particle indices in their own contiguous array so the head
   walk reads sequentially; compute each particle's `npeaks` once. The
   frontier loop keeps its order and its single `ran3()` draw per iteration,
   the swap-remove order of `sync_frontier_class` and the unstable `hpsort`.
   A tester keeps the present routine as the reference and compares the
   assignments on synthetic sparse tables, ties included, at 1 and 11
   threads.
5. **The compact dense 2D record.** As 3D: 8 bytes per entry (distance and
   in-plane index; the class is the row index), plus a list of refined
   entries (particle, class, shift) for the entries whose shift was refined
   (at most `nhood_sz` per particle); membership in the list stands for
   `has_sh`. A format word and a class-count check in the header. The master
   holds the rows and the refined list in compressed-row form and rebuilds a
   full candidate through an accessor; the old 2D dense reader goes (part
   files live for one iteration; no older format is read). The 2D tester
   writes the old format and sets `loc_tab` directly: change it with the
   format. Teach the lossless part-file decoder the new layout and compare
   every field both layouts hold, bit for bit.
6. **The frontier loop's sort** (only if item 1 shows the loop dominates
   after items 2-4): a top-k selection instead of sorting about 150 values in
   `sample_likelihood_index`, falling back to `hpsort` whenever equal values
   reach the top k+1. A tester compares it with the present sampler on
   random and tied inputs, draw for draw.

### Phase 2: 2D look and feel (outputs unchanged; log text and interface text change)

7. **Log lines as in 3D.** One `>>> CONVERGED: yes|no` per iteration with the
    final decision, after `minits`, `maxits` and coverage (the decision logic
    unchanged); `>>> % PARTICLES USED FOR UPDATING THE MODEL`; `>>> MATCHING
    LOW-PASS LIMIT (THIS ITERATION)` from the particles' `lp` field; the
    meaningless `PROB NEIGH MODE` line and the spread of a constant low-pass
    gone; the stage line `>>> STAGE n WITH LP x A` (format fixed) and the
    gold-standard stage line in the same style; `>>> SOLVE2D ENTERING REFINE2D
    STAGE n` and `>>> SOLVE2D RETURNED FROM REFINE2D STAGE n` with a flush;
    `prob_align2D` and the 2D assignment one line per iteration; the stage
    table without the `MAXCLSPOP` column, and the dead `max_cls_pop` deleted.
    The stats file's keys stay. `check_conv2D` has two callers, the 2D matcher
    (shared memory, and the stream's 2D pool through it) and the refine2D
    strategy (distributed); each prints the final verdict after its own
    `minits`/`maxits`/coverage decision, so the stream pool keeps one
    convergence line per iteration without any edit to stream code.
8. **The solve2D interface.** Correct help texts: `center` `{yes}`;
    `lpstart`, `lpstop`, `cenlp` described as derived from the mask diameter
    (with the rule); `nsample` Standard, showing its real default; `nstages`
    exposed with the 2D meaning (it re-spaces the low-pass ladder and moves
    the gold-standard stage); `refine=snhc_smpl` explained (it differs from
    the default only in the last stage). Delete the unused `refine2D`
    `ui_program`. Teach `scripts/default_audit.py` to follow `exec_solve2D`
    into `exec_solve2D_workflow` if that is a small change; otherwise record
    the generator gap for the maintainer. Check the generated defaults and
    the UI JSON before and after; only the intended texts change.
9. **Class-average names through helpers.** A `simple_refine2D_fnames.f90`
    beside `simple_refine3D_fnames.f90` (same style: typed `string`
    functions on the `simple_defs_fname` constants) for the merged, even,
    odd, ranked and JPEG class averages of an iteration, the extension an
    argument; every non-stream site converted, including the parse of the
    iteration in `make_cavgs`; the dead `refine3D_iter_refs_fname` deleted.
    Every file name stays byte-identical (identity compares names). Stream
    sites listed in the report.
10. **Shared sigma2 helpers.** `prepare_canonical_sigma_update` and
    `record_sigma2_bench_context` in one module used by both strategies (and
    by the other variants where the body is the same); `ensure_canonical_sigma2`
    shared with explicit arguments for what differs, if that stays readable;
    otherwise leave it and say why.
11. **`solve2D_policy.md`** gains Defaults and Stage-controller sections like
    `solve3D_policy.md`; its sigma2 line describes the one canonical store and
    its generation transaction; the "final fill-in pass" lines match the
    code (`fillin=yes` on the staged calls).

### Phase 3: the 3D table, the shift-search kernels and the NU filter

12. **The 3D master table by reference.** Store the distances so that each
   reference's distances over the particles are contiguous (for example
   `dists(nptcls,nrefs)` and `inpls(nptcls,nrefs)`), read consecutive local
   particles of a part in one statement and scatter them through the
   global index, and place refined entries at their final offsets during the
   read (counts per particle from the parallel read, a prefix sum, then
   placement in a fixed order; the testers compare the refined arrays element
   by element). Restructure `assign_state_labels` and the check loop with
   references outside and per-particle best-so-far arrays, keeping the strict
   `<` in ascending reference order. Delete the dead `score_min`,
   `score_spread` and `score_mode` of `prepare_ref_score_vectors`. Every index
   into the table changes (`simple_eul_prob_tab.f90` and its tester); the
   part-file format does not. Stress test of the read.
13. **The shift matrix fused into the four alignment kernels**
   (`gen_corr_cc_for_rot_8_2`, `gen_euclid_for_rot_8_2`,
   `gen_corr_cc_grad_for_rot_8`, `gen_euclid_residual_grad`): a per-thread
   vector for the current shell, filled with the first shell and the per-shell
   step by the same statements as `gen_shmat4aln_8`, used where the kernel
   reads `shmat_8(rp,k)` and advanced after each shell. `gen_shmat4aln_8`
   stays for `calc_frc` and `gen_sigma_contrib`. A tester compares old and new
   kernels bit for bit (value and gradient, `irot` at and beyond `pftsz`,
   both objectives, several shifts); if a single bit differs (the compiler's
   fused multiply-adds may), restore the item and record it. Delete the dead
   store `denom = self%wsqsums_ptcls(i)` in the correlation gradient.
14. **The NU label optimisation.** Sub-timers inside `optimize_nu_cutoff_finds`
   first (level smoothing against the Potts sweeps). Then: reuse the setup's
   smoothed candidates and its packed norms at each level (42 smoothings to
   28); skip all-zero lines in the tent smoothing, writing their zeros
   (detect them from the data: `simple_image_filt.f90` uses the same routine);
   block the strided y and z passes, each column keeping its own operation
   sequence; packed voxel lists per colour in the Potts sweeps, the neighbour
   sum kept in ascending order with its `/degree`. Bitwise tester checks of
   the tent smoothing and of the label field against the present code, in
   `simple_nu_filter_tester.f90`. NU3D and the cases with the filter judge the
   item.

### Phase 4 (only if at least two and a half hours remain): 2D cropped particles shared

15. **The 2D table on the cropped image.** One shared function decides
    eligibility for the 2D table and the 2D matcher (not cached, not the
    stream's 2D pool, `box_crop < box`, `ptcl2D`). When eligible, `prob_tab2D`
    makes the in-memory cache entry of each particle (noise-normalized at the
    full box, Fourier-cropped, unshifted, back in real space) and prepares
    alignment from it through `prepimg4align_cached`, without needing
    `ptcl_imgs`; the record is kept (or put into the buffer) before alignment
    is prepared from it, because `fft_clip_shift` transforms its input in
    place. In shared memory the table opens the particle buffer for records
    (`open_record_buffer`, `put_record`, `seal_record_buffer`) instead of
    filling it with raw images, and the matcher takes a batch's records
    (`take_records`, then a parallel copy into `ptcl_imgs`) and skips the read
    and the full-box preparation; otherwise it works as today. In distributed
    runs a table worker then does a little more per particle (a `box_crop`
    inverse and forward transform); the report says so. Identity: E2 and E2S
    with `cache=no` on the item's build must be byte-identical to E2C and E2SC
    (`cache=yes`) on the base build, apart from the files Phase 0 declared for
    this comparison by reason (the cache files, and the files that record the
    command line, which says `cache=yes` in one and `cache=no` in the other).
    Any other difference: the item is restored. Stress test of the parallel
    copy.

## The identity test

As in the runs `speed_mac` and `speed2d_mac` (`identical_speed_run.md`,
section "The identity test"; `speed2D_run.md`): one-thread runs with a fixed
seed, each case in a fixed directory, outputs moved aside after each run,
every regular file compared except exclusions declared in Phase 0 with
reasons, project files compared through `print_project_field` dumps (ptcl2D,
ptcl3D, cls2D, cls3D, out), table part files compared through the lossless
decoder (accepted by the maintainer; teach it each new layout). The fixture's
checksums are checked after every case. The cases run in the CTest
environment with `OMP_NUM_THREADS=1`; `SIMPLE_SEED` is set for `simple_exec`
cases. The harness scripts of `~/speed2d_mac/bin` and `~/speed_mac/bin`
(`run_case.sh`, `compare.py`, `probtab_decode.py`, `item_src.sh`,
`snap_item.sh`) may be copied after checking them, if those directories still
exist; otherwise write them again.

**Cases.** All made again on this run's base build.

| case | what it covers |
| --- | --- |
| SW1, SW6 | the simulated workflows (1jxy with its two-part step, 6vxx), including their `solve2D` steps |
| MSI, CCI, PCI | the three-state `refine3D` (dense then neighbourhood table), the correlation objective with fractional update and trailing reconstruction, the PCG backend (from the SW1 fixture, as before) |
| NU3D | `nu_filt3D` and `postprocess_nu` on the fixture's final half maps (the filter's standalone programs) |
| E2 | `solve2D` with two parts at one thread on about 1,500 particles simulated from the fixture's 1jxy truth map (CTF, random poses and shifts), about 20 classes, all six stages: the `snhc_smpl`, sparse `prob_snhc` and dense `prob` table paths, class-average assembly, the distributed master |
| E2S | E2's command line in shared memory (the shared-memory iteration, the particle buffer) |
| R2 | `refine2D refine=prob` for two iterations from E2's final project, two parts (the dense 2D table, its reads and its assignment) |
| E2C, E2SC | E2 and E2S with `cache=yes` on the base build only: the reference for Phase 4 |
| AD1, CR1 | `solve3D_addon` and `cont_refine3D_1jxy` at one thread, once at the start and once at the end of the run, compared by their check values |

E2 and E2S must crop (`box_crop < box`), so that the cropped class sums and
the records are exercised; check it in Phase 0 from the logs or reports and
size the particle set so they crop and each 2D case takes at most about ten
minutes at one thread. Check that the 2D cases reach every stage and table
path they claim. Check that MSI, CCI and PCI report the cropped records in
their refine3D reports (the line `match3D particle buffer` padded, then `:
cropped records`; grep `particle buffer +: cropped records`). Check that no
compared file holds log text: items 1 and 7 change log lines, and the
exclusions are frozen after Phase 0.

**The 3D record buffer of `af13bc828`.** That commit made the 3D
reconstruction take its observations from records kept by the search, and
claims byte-identical results; it was never run. In Phase 0 build its parent
`aea5b1846` as well (`git --no-optional-locks archive aea5b1846 | tar -x -C
<dir>` in the scratch directory) and run MSI, CCI and PCI on it once: they
must equal the base runs byte for byte (these cases run no 2D step, so the
commit's 2D change cannot reach them). A difference is a defect of that
commit: fix it as the first item of Phase 0, or, if it cannot be fixed
identically, stop with the evidence.

**When to run what.** After every item: the item's Debug build and the
fast-gate suites, the identity cases the item can reach (SW1, E2, E2S and R2
for 2D items; SW1, MSI, CCI, PCI and NU3D for 3D and filter items; all of them
for shared code), and the stress tests of any parallel region it touches. At
the end of every phase: SW6, the full `./compile_debug.sh` on an idle Mac, and
one multi-threaded run of E2S at 11 threads and of MSI with two parts of 5
threads (they must run to the end with finite results; not byte-compared).
AD1 and CR1 at the end of Phase 3, or at the start of Phase 5 if Phase 3 did
not reach its end (section "The clock").

## Common rules

As in `identical_speed_run.md`, section "Common rules", with these changes:

- The base is master `af13bc828` or later, whatever `~/src/SIMPLE` holds at
  launch; the worktree is `~/src/SIMPLE-wt-speed2d3d`.
- Scratch is `~/speed2d3d_mac`; the previous runs' scratch directories may be
  read, never written.
- No new test suites or registrations: new checks go into existing tester
  modules of registered suites (`simple_eul_prob_tab2D_tester.f90` and
  `simple_eul_prob_tab_tester.f90` for the tables, `simple_pftc_inplane_tester.f90`
  for the kernels, all in `unit_pftc_align2D3D`; `simple_nu_filter_tester.f90`
  in `unit_reconstruction` for the filter and the tent smoothing). New source
  files under `src/` are picked up by the existing globs; no CMake edits. The
  globs have no `CONFIGURE_DEPENDS`: after adding a file, configure again
  (`cmake -S . -B build`, and for each Release build) before building.
- Fix as you go; an item that cannot be made identical is restored and
  recorded; per-item snapshots and diffs; nothing committed.

## Phases

**Phase 0: base and cases.** The Debug build with the fast gate, the base
Release build of the worktree as it is, and the parent build of
`aea5b1846`. The fixture (SW1 once, in a fixed directory), then the E2
particle set, then every case twice on the base build (AD1 and CR1 once each,
as a reference for the end of the run; E2C and E2SC twice, for Phase 4), and
MSI, CCI and PCI once on the parent build. Wall time and peak memory
recorded. The exclusion list from the pairs, and, declared by reason, the
files the Phase 4 comparison of a `cache=no` run with a `cache=yes` run may
skip (the cache files; the files that record the command line). Exit: every case deterministic at one thread
under the list; the 2D cases crop and reach their paths; the parent and base
agree on MSI, CCI and PCI.

**Phase 1** (items 1-6), **Phase 2** (items 7-11), **Phase 3** (items 12-14).
Exit of each: every kept item identical on its cases and passing its stress
tests, its diff written; items not kept or not reached recorded with the
cause; the phase-end runs done.

**Phase 4** (item 15), only with at least two and a half hours left;
otherwise its status is `DONE: not started, for time`.

**Phase 5: documents and close.** `doc/policies/2D/solve2D_policy.md` (the 2D
table formats, reads and assignment; the look and feel; Phase 4 if done),
`doc/policies/3D/refine3D_policy.md` and
`doc/policies/NU/nonuniform_filtering_policy.md` where they describe what
changed, `particle_read_path_policy.md` and `particle_cache_policy.md` if
Phase 4 ran; `speed_architecture_proposals.md`: what is now done (by this run
and by `speed2d_mac`, whose items it still lists as open), and the two
corrections of "What the reading found" (the 3D table is already
particle-major; the shift-independent term is not identical); this note's
progress; the report beside it (`speed2d3d_run_report.md`): what changed per
item, the identity and stress evidence, the time indications from item 1 and
the phase-end runs before and after, what was not kept or not reached and
why, proposed skill text, and a commit message. Then `full.diff` in the run
directory.

## File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/completed/speed2d3d_run.md` | this plan, its progress, the report beside it | each |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90`, `simple_eul_prob_tab_utils.f90`, `simple_eul_prob_tab2D_tester.f90` | items 2-6, 7 | 1 2 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | items 1, 3, 7, 15 | 1 2 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `simple_eul_prob_tab_tester.f90` | item 12 | 3 |
| `src/main/pftc/simple_polarft_corr.f90`, `_calc.f90`, `_geom.f90`, `simple_pftc_inplane_tester.f90` | item 13 | 3 |
| `src/main/nu_filt/simple_nu_filter_bank.f90`, `_state.f90`, `_potts.f90`, `_vars.f90`, `_tester.f90`, `src/utils/filter/simple_tent_smooth.f90` | item 14 | 3 |
| `src/main/simple_convergence.f90` (`check_conv2D`), `src/main/commanders/simple/simple_commanders_solve2D.f90`, `src/main/solve/simple_solve2D_controller.f90`, `src/main/strategies/search/simple_strategy2D_matcher.f90` | items 7, 9 | 2 |
| `src/main/ui/simple/simple_ui_refine2D.f90`, `scripts/default_audit.py` | item 8 | 2 |
| `src/defs/simple_refine2D_fnames.f90` (new), `src/defs/simple_refine3D_fnames.f90`, the class-average name sites (refine2D strategy, 2D matcher, `simple_commanders_mkcavgs.f90`, `simple_starproject.f90`, `single_commanders_nano2D.f90`) | item 9 | 2 |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90`, a shared sigma2 helper module | items 7, 9, 10 | 2 |
| `src/main/strategies/search/simple_matcher_ptcl_batch.f90`, `simple_matcher_ptcl_io.f90`, `simple_strategy2D_matcher.f90`, `simple_matcher_2Dprep.f90`, `probabilistic/simple_eul_prob_tab2D.f90` | item 15 | 4 |
| `doc/policies/2D/solve2D_policy.md`, `doc/policies/3D/refine3D_policy.md`, `doc/policies/NU/nonuniform_filtering_policy.md`, `doc/policies/particle_read_path_policy.md`, `doc/policies/2D/particle_cache_policy.md`, `doc/implementation_notes/planned/speed_architecture_proposals.md` | documents | 2 5 |

## Progress

### Phase 0, 2026-10-19: base and cases

The base is master `af13bc828`; scratch is `~/speed2d3d_mac`, the run
directory `~/agent_runs/runs/speed2d3d_mac`. No source file changed.

- **Builds.** `./compile_debug.sh` passed its fast gate, 15 of 15 suites
  in 14.8 s. The base Release build and the parent build (`aea5b1846`) are
  in `~/speed2d3d_mac/builds` (gfortran 16 from Homebrew, Release,
  `BUILD_TESTS=ON`, OpenMP, `USE_ARCHOPT=ON`).
- **The 3D cases crop.** The fixture's box is 192, and `refine3D` crops
  only when `box_crop` is given, so MSI, CCI and PCI take `box_crop=128`.
  Their refine3D reports then say `cropped records`.
- **The 3D record path of `af13bc828` is identical.** MSI, CCI and PCI on
  the parent build equal the base byte for byte (216, 50 and 53 files).
- **The 2D cases.** 1,500 particles simulated from the fixture's 1jxy map
  (box 192, 1.3 A). E2 and E2S crop to 88 pixels and reach all six stages:
  `snhc_smpl`, the sparse `prob_snhc` table, the dense `prob` table, class
  averages assembled by the master. R2 runs two dense iterations on a copy
  of a fixed E2 project. E2C and E2SC are E2 and E2S with `cache=yes`; E2CN
  and E2SCN run E2 and E2S with `cache=no` in the directories of E2C and
  E2SC, because `metadata.json` and the project dumps hold the case path.
- **Determinism.** Every case ran twice with no difference under the
  declared exclusions (`exclusions.tsv` in the run directory): SW1 614
  files, SW6 546, MSI 216, CCI 50, PCI 53, NU3D 14, E2 140, E2S 137, E2C 140,
  E2SC 137, R2 30. AD1 (46 check values) and CR1 (60) pass. E2S at 11
  threads and MSI with two parts of 5 threads ran to the end.
- **Cost** (one thread, wall time, peak memory): SW1 11 min, 4.4 GB; SW6
  13 min, 3.8 GB; MSI 2 min, 1.6 GB; CCI 12 s; PCI 42 s, 1.7 GB; NU3D 1 min;
  E2 1.5 min, E2S 2.3 min, R2 17 s; AD1 18 min, 2.1 GB; CR1 33 min, 10.8 GB.

### Phase 1, 2026-10-19: the 2D tables and the master's work between them

Five items kept, each built from its own snapshot and byte-identical to the
base on SW1, E2, E2S and R2 (the part files R2 keeps compared through the
lossless decoder). They were done in the order 1, 2, 5, 3, 4, so that the
parallel read was written once, for the final formats.

- **Timers.** The `PROB_ALIGN2D` report has a section that splits the reads
  (the parts' heads, their candidates) and the assignment (class lists:
  counts, fill, sort; the frontier loop; the leftover particles).
- **One pass over each sparse part.** The 2D part file of the candidate
  stream has a format word (-4) and the candidate count of every particle
  after the particle indices. The master reads every part's head, sizes
  its store from the counts, and reads every candidate once.
- **The compact dense record.** The dense part file (format -5) holds one
  row of (distance, in-plane index) per particle and class, 8 bytes per
  entry, then the refined entries (particle, class, shift). The master
  keeps the rows and the refined entries in compressed-row form and
  rebuilds a full candidate when it assigns. R2's part file shrinks from
  432 to 155 kB.
- **Parallel reads.** Both passes run in parallel over the parts. Each
  part's size is checked against its header, and a particle in no part or
  in several stops the run. Errors are raised after the loop.
- **The sparse assignment's class lists.** Each class's candidates, with
  their particles and distances, are gathered during the fill, sorted in
  parallel over the classes with the same unstable `hpsort`, and walked
  sequentially. Peak counts are computed once. The frontier loop, its one
  draw per iteration and its tie order are unchanged.
- **Tests.** The 2D tester writes its part files through the table's own
  writer, checks that a refined shift survives the dense layout, reads 7
  parts at one thread and 20 times at every processor (both layouts), and
  compares the sparse assignment with the algorithm it replaced on tables
  with ties, at four draw seeds and at one thread against every processor.
- **Not done.** The top-k sampler of the frontier loop: on E2 every phase
  of the read and assignment takes 1 to 3 ms of a 2 s iteration, so the
  frontier loop does not dominate. The timers settle it on a large case.

### Phase 2, 2026-10-19: 2D look and feel

Five items kept. The code items were each built from their own snapshot and
are byte-identical to the base on SW1, E2, E2S and R2 (the sigma2 helpers
also on MSI, CCI, PCI and NU3D); only log and interface text changed.

- **Log lines as in 3D.** `check_conv2D` tests the criteria and no longer
  prints a verdict; the 2D matcher and the distributed refine2D strategy
  each print one `>>> CONVERGED: yes|no` after their own `minits`, `maxits`
  and coverage decision. The stream's 2D pool keeps one line per iteration
  through the matcher, with no edit to stream code. The low-pass line reads
  the particles' `lp` field, `PROB NEIGH MODE` and the spread of a constant
  low-pass are gone, and the stage lines, the entering and returning lines,
  and one line each for `prob_align2D` and the 2D assignment match 3D. The
  stage table lost `MAXCLSPOP`, and the dead `max_cls_pop` is deleted. The
  stats file's keys are unchanged. E2's log went from 872 to 760 lines.
- **The solve2D interface.** Correct help texts for `center`, `refine`,
  `lpstart`, `lpstop` and `cenlp` (with the mask-diameter rules); `nsample`
  is Standard and shows 200,000; `nstages` is exposed with its 2D meaning;
  the unused `refine2D` program variable is deleted. `default_audit.py`
  follows `exec_solve2D` into `exec_solve2D_workflow`, so the interface now
  shows the defaults the code applies (`autoscale=yes`, `cls_init=rand`).
  Only solve2D's entry of the UI JSON changed. The numeric defaults shown
  for `lpstop` and `cenlp` are still parameter baselines: the interface has
  no way to show a derived default (a question for the maintainer).
- **Class-average names through helpers.** `simple_refine2D_fnames.f90`
  builds the merged, half and ranked class-average names of an iteration
  and parses the iteration back; every non-stream site uses it. The stream
  still builds five names inline (`simple_stream_pool2D.f90`,
  `simple_stream_refine2D_utils.f90`).
- **Shared sigma2 helpers.** `simple_strategy_sigma2.f90` holds the store
  preparation and the bench counts both strategies copied. The strategies'
  `ensure_canonical_sigma2` stay apart: their shared core is already one
  call, and the rest is what differs.
- **`solve2D_policy.md`** has Defaults and Stage-controller subsections,
  describes the canonical sigma2 store, and states the fill-in as the code
  does it (`fillin=yes` on the sampled last stages).

### Phase 3, 2026-10-19: the 3D table, the shift-search kernels and the NU filter

Two of the three items kept; the kernel fusion is restored.

- **The 3D master table by reference.** The master holds the dense table as
  `dists(nptcls,nrefs)` and `inpls(nptcls,nrefs)`, so each reference's
  distances are contiguous. A part's consecutive rows are read in one
  statement (up to 8 MiB) and scattered to the columns; the read counts
  each particle's refined entries, and every part places its entries at
  their final offsets in parallel. The invalid/flat check and the state
  labels run references outside, in blocks of 256 particles, keeping the
  strict `<` in ascending reference order; the dead `score_min`,
  `score_spread` and `score_mode` are gone. Byte-identical on SW1, MSI,
  CCI, PCI and NU3D; the tester's 24 parallel reads at every processor
  match the read at one thread.
- **The shift matrix fused into the alignment kernels: not kept.** The four
  kernels advanced one shell vector per thread with the same statements as
  `gen_shmat4aln_8`, yet the results differ (one SW1 particle's shift in
  the 7th digit; E2, E2S and R2 broadly): the compiler's fused
  multiply-adds differ between the two contexts. Restored, as the plan
  says.
- **The NU label optimisation.** Sub-timers split the optimisation into the
  level smoothings and the Potts sweeps (`VOLASSEMBLE` report). Each level
  takes its own candidate and its smoothing norm from the setup instead of
  smoothing them again (42 smoothings to 28). The tent smoothing writes
  all-zero rows, planes and slabs as zeros and filters the Y and Z columns
  a plane at a time, each column keeping its operation sequence; a tester
  compares it bit for bit with the column-by-column version it replaced,
  at one thread and 20 times at every processor. The Potts sweeps walk a
  packed voxel list per colour; a tester repeats the optimisation 20 times
  at every processor and finds the same label field. Judged together on
  SW1, MSI, CCI, PCI and NU3D (the parts' own snapshots held the kernel
  change).

### Phase 4, 2026-10-19: 2D cropped particles shared

The item is kept.

- **The 2D table on the cropped image.** One function, `crop_particles2D`,
  decides for the 2D table and the 2D matcher (uncached, outside the
  stream's 2D pool, `box_crop < box`, `ptcl2D`). Then the table makes each
  particle's in-memory cache entry (noise-normalized at the full box,
  Fourier-cropped, unshifted, back in real space) and prepares alignment
  from it as a cached particle, as the matcher already did.
- **Records in shared memory.** The table opens the particle buffer for
  records instead of filling it with raw images and puts each cropped
  particle before alignment transforms its copy. The matcher of the same
  iteration copies a batch's records in parallel and neither reads nor
  prepares it at the full box: in E2S's 14 table iterations it read no
  particle. A distributed table worker crops each particle itself.
- **Identity.** `solve2D` with `cache=no` on the item build is
  byte-identical to `cache=yes` on the base build, distributed (140 files)
  and in shared memory (137 files). Uncached results therefore change by
  rounding (accepted with the package). E2S ran 20 times at 11 threads,
  each to the end with finite class averages.

### Phase 5, 2026-10-19: documents and close

- `solve2D_policy.md`: the 2D part-file layouts, the parallel reads, the
  sparse assignment, the cropped table and its records.
  `particle_read_path_policy.md` and `particle_cache_policy.md`: the records
  of the 2D table, and uncached 2D equal to cached. `refine3D_policy.md`: the
  table by reference and its read. `nonuniform_filtering_policy.md`: the
  reused smoothings, the tent smoothing, the packed Potts sweeps, the
  sub-timers.
- `speed_architecture_proposals.md`: a status section (what this run and
  `speed2d_mac` did) and the two corrections to section 1.3.
- The report: `speed2d3d_run_report.md`.
