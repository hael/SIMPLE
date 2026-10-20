# Identical speed changes: the Mac run `speed_mac`

Implementation note, 2026-10-17. **Carried out 2026-10-17** (progress
below; report in `identical_speed_run_report.md`). The
execution plan of the unattended run `speed_mac` on the maintainer's Mac,
made while the Oracle Linux test machine runs `test_truth`. It carries out
the parts of the speed week that leave every result byte-identical: item 3
of `multistate_speed_run_completed.md` (the table reads), days 3 and 4 of
`single_state_speed_week_completed.md`, and the identical items of that plan's section
"Multi-state". The design is in those two notes; this note holds what the
code check of 2026-10-17 changed, the identity test, the phases and the
progress.

## Rulings of 2026-10-17

- The run is a headless Claude Code session on the Mac, in a worktree of
  master, launched by the maintainer.
- **Identical changes only.** Every item must leave the outputs of the
  identity cases byte-identical. Changes that alter results (single-precision
  weights in the score kernels, a cap on refined references, the worker fold
  of day 5) are not in this run.
- **Identity now, timings later.** The run proves identity. Wall-time gains
  are measured afterwards on the test machine (once `test_truth` is done) and
  on the maintainer's PfCRT and msp1 runs. The run may time single routines
  where an item's gain is in doubt (below).
- The run may use most of the Mac (11 cores).

## What the code check of 2026-10-17 changed

Five read-only checks of the code at `25d256886` (one per area) changed the
plan in these places. Line numbers are approximate.

**Table reads** (`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`,
`simple_commanders_prob.f90`).

- The master fills its dense table twice before it reads a part: once by
  `allocate`, because `prob_candidate` has default values
  (`simple_eul_prob_tab_utils.f90`, about line 20), and once in
  `initialize_storage`. Both passes sit inside the reported read time.
- The initialisation cannot simply be dropped. The neighbourhood table uses
  the default values for the placeholder entry of a particle without
  candidates (`read_tabs_to_glob`, about line 1283), and the table tester
  sets only `dist` and relies on `initialize_storage` for the rest (tester,
  about lines 95-100 and 135-140). The dense table therefore needs its own
  entry type without default values, which the compact record (below)
  provides.
- An OpenMP loop over the part files is safe: each part writes only its own
  particles' columns, and every buffer and unit is local. One exception: the
  read-modify-write of `self%seed_nrots` (about lines 932-933) races; collect
  the per-part values and check them after the loop.
- The read does not check that every particle of the master appears in
  exactly one part. Today a missing particle keeps the initial values and an
  overlap lets the last part win silently. The parallel read adds that check.
- In shared memory there is one part file, so the parallel read gains
  nothing there.
- `read_tabs_to_glob` does open every part twice: the first pass reads the
  particle indices of every chunk to count the candidates per particle (about
  4 of every 32 bytes), the second places the candidates. The order matters:
  the best evaluated reference is chosen with a strict `<`, and the
  assignment graph lists each reference's particles in file order. Reading
  every part into memory first would hold the candidates twice on the master,
  so the one-pass read is done only without that copy: the part writes its
  per-particle candidate counts after the particle indices, the first pass
  reads only those, and the second pass is unchanged. The gain is modest; if
  this is not simpler than today's code, the item is dropped and the reason
  recorded.

**The compact record for the dense table.**

- An unrefined entry is always (reference, in-plane index, distance, 0, 0,
  no shift). A refined entry replaces distance and in-plane index and sets
  the shift, only when the search returns `irot > 0`.
- `ref_assign` reads only `dist`; `assign_candidate` reads all fields of the
  assigned entry; `simple_strategy3D_prob.f90` reads only the assignment map.
  Nothing reads the shift of an unrefined entry.
- So a row of (distance, in-plane index) per particle and reference, 8 bytes,
  plus a short list of (reference, shift) per particle for the refined ones,
  reproduces the same assignment and the same values handed on. File entries
  go from 28 to 8 bytes and the master's from 24 to 8.
- The neighbourhood table shares the writer (`begin_write`,
  `flush_candidate_buffers`, `write_tab`) and the empty-part writer
  (`write_empty_prob_tab`). The dense table gets its own writer, reader and
  empty-part writer, with a format word in the header; the neighbourhood and
  2D tables keep theirs.
- Side gain: `append_or_replace` (utils, about line 113) scans a particle's
  entries backwards for every refined reference; with rows indexed by
  reference it is one store.

**The shift refinement** (`src/main/pftc/`).

- Every L-BFGS-B evaluation goes through `gen_corr_grad_for_rot_8`
  (`simple_polarft_corr.f90`, about lines 881-900), which copies the
  reference to double precision; the kernel then rebuilds the shift matrix
  (`gen_shmat4aln_8`).
- The "copied twice" path of the week plan is dead:
  `gen_corr_grad_only_for_rot_8` is called only from `grad_shsrch_gcostfun`,
  which L-BFGS-B never reaches because the shift search sets the combined
  cost-and-gradient function. It is deleted. `gcostfun_8` stays set, since
  `simple_opt_lbfgsb.f90` checks that one gradient function is associated.
  `shift_ref` (`simple_polarft_geom.f90`) has no caller and is deleted.
- The reference does not change during one minimisation, and the conversion
  to double precision is exact, so one copy per minimisation gives the same
  bits. The copy is about a tenth of one evaluation: 3-8% of the refinement.
  It goes into a dedicated per-thread buffer in the heap variables
  (`simple_polarft_core.f90`), loaded at the start of `minimize`; the kernels
  themselves stay as they are. `pft_ref_8` cannot be reused (`calc_frc` and
  `gen_sigma_contrib` overwrite it), and a cache keyed on the reference index
  would serve stale data when references are replaced under the same index.
- **New: a repeated objective.** After L-BFGS-B the alternation loop of
  `grad_shsrch_minimize` (`simple_pftc_shsrch_grad.f90`, about line 397)
  evaluates the objective over all rotations at the final shift. The
  L-BFGS-B callback evaluates the same objective at every accepted point
  (`simple_opt_lbfgsb.f90`, about lines 143-150), so when the minimiser ends
  on an accepted point the vector is recomputed bit for bit. It is not always
  so: in the multi-state table the start rotation comes from the
  particle-frame grid, not from a vector at the start shift, and if L-BFGS-B
  stops at its start point the loop can take a second pass. So: keep each
  all-rotation vector computed by `gen_objfun_vals` (before the loop, in the
  callback, in the loop) with its shift, and reuse it only when the shift is
  bitwise equal; the loop stays. This is worth more than the copy.
- **The shift-matrix range:** `gen_shmat4aln` fills the double-precision
  matrix up to `kfromto(2)` and then converts it up to `interpklim`; the
  shells above `kfromto(2)` hold leftovers. Every consumer loops over
  `kfromto` only. Limit that one conversion to the searched shells; leave
  `gen_shmat` (used by `shift_ptcl`, all shells) alone.
- **New: allocations in the score loop.** `gen_likelihood_val_bank` allocates
  `scores(nrots)` per reference (corr, about line 576) and
  `sample_likelihood_dist` allocates its `cdf` (utils, about line 490); the
  distance function is an internal procedure using host variables (corr,
  about lines 593-601), which gfortran may implement with a trampoline per
  call (`-Wtrampolines` tells). Per-thread buffers and a sampler working on
  an array of distances, evaluated in the same order, give the same results.
- What must not change for identity: the order of the floating-point
  operations in the shift matrix (product, then sine and cosine of the first
  shell and of the step, then repeated multiplication) and in the kernels
  (the split of the angle loop, the running sums per shell, the weighted sum
  and the normalisation). Release builds use `-O3` (with `-mtune=generic` on
  the Mac, `-march=native` on the test machine) and fuse multiply-adds, so
  any rewrite of those expressions can change bits.

**The nonuniform filter** (`src/main/nu_filt/`, `simple_nu_state_filter.f90`).

- Every candidate low-pass (up to 8) is filtered, written as an MRC pair,
  read back for the objective and read again for the composition. The write
  and read are lossless (float32, no rescaling), and nothing else reads the
  cache files.
- The objective already runs one candidate at a time; the optimiser needs
  only the packed distance maps, which are held anyway. So the objective can
  run right after each candidate is filtered, with no volumes kept.
- The composition already has a path that filters again from the input pair
  (`compose_nu_labels_from_volume`, `simple_nu_filter_apply.f90`, about lines
  112-158): the same copy, taper, forward FFT, Butterworth multiply and
  inverse FFT, with plans made by `FFTW_ESTIMATE`, so the voxels are
  bit-identical. The input pair is killed just before the composition
  (`simple_nu_state_filter.f90`, about lines 137-147); moving the kills below
  it is enough. The composition's peak (about 6 volumes) stays below the
  setup's present peak (about 8).
- `postprocess_nu` writes the caches and never composes; `resolest` composes
  from them with the input pair still alive. Both follow.
- There are no timers inside the filter: the volume-assembly report has one
  line for all of it. The optimiser's tent smoothings (35 at 8 candidates)
  may cost more than the disk. Filtering again in the composition trades
  reads from the page cache for two inverse FFTs per candidate, which may not
  pay. Timers come first, and the composition change is kept only if a
  timing of the filter on a realistic pair says it is not slower.

**The states' volume assembly in parallel is not an identical change**, and
it leaves the run.

- Several sums in a state's assembly depend on the thread count: the radial
  FSC (a double-precision reduction), the noise term of the ML
  regularisation, the smoothing weight of the filter's label update (a
  single-precision reduction that can flip near-tie labels), the phase
  randomisation (per-thread random streams), and FFTW plans made for the
  thread count. States run concurrently with fewer threads each would
  therefore differ from the serial loop, beyond rounding in the labels.
- The filter is module state, its cache names carry no state index, and
  `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
  forbids holding arrays of state volumes in assembly.
- Instead, three identical items of the same code: per-state timers in the
  volume-assembly report; deleting two dead volume reads per state
  (`restore_state_from_parts`, `simple_commanders_rec_distr.f90`, about lines
  111-112, and the stale comment at 44-45); and filling the filter's source
  halves from memory instead of reading back the files just written (about
  lines 460-479; copy the logical region only, since an image copy carries
  the FFT padding that a read leaves at zero).

**The reprojection model and the project reads** (the distributed
`refine3D` strategy).

- `materialize_reprojection_model` (`simple_refine3D_strategy.f90`, about
  lines 525-553) builds a builder every iteration: a full project read,
  the toolboxes, the FSC read, the projection, and a write of the particle
  segment. The toolboxes are cheap except the symmetry distance table in
  neighbourhood mode; the project read is the cost.
- Crop, directions and masks are fixed for one `refine3D` run (a stage that
  changes the crop starts a new one); the project changes every iteration.
  A persistent builder therefore keeps its toolboxes and refreshes its
  project every iteration, keeps the per-iteration `set_projs` check, and
  never holds a second full project between iterations. It is an allocatable
  component of the strategy, released in `distr_cleanup` (compile-time
  policy).
- `prob_align`, volume assembly, the post-processing and `commit_sigma2` run
  in the master process. Per distributed iteration the master reads the full
  project five times (the reprojection model, the strategy, `prob_align`,
  volume assembly's `initialize_context`, and `commit_sigma2`, which the week
  plan missed) and the particle segment twice (volume assembly's trailing
  reread, and `refresh_matching_lp_from_project`, which uses one record).
- Two reads can go at once: volume assembly's trailing reread (nothing
  writes in between), and the one-record read (`read_segment` with
  `fromto=[fromp,fromp]`).
- Handing the in-memory project to `prob_align`, volume assembly and
  `commit_sigma2` needs a new entry, since commanders take only a command
  line; the PCG assembly from the master's `(params, build)` is the
  precedent. They must work on copies, not in place (the master never takes
  volume assembly's resolution fields, and writes the whole project at the
  end), and `commit_sigma2`'s copy must reproduce what `build_spproj` sets
  for the project's `cwd`. Whether copying a 750,000-particle project is much
  faster than reading it is unmeasured; the run measures it first.

## The identity test

An identical change is judged only by identity: the same case on the build
before and the build after, at one thread with a fixed seed, outputs
compared byte for byte. Multi-threaded runs are never used to judge an
identical change; their draws vary from run to run (the run `ms_speed`
stopped on such checks).

**The cases.** All run with `OMP_NUM_THREADS=1`, in the environment of the
CTest registration (`production/CMakeLists.txt`, `simple_add_test`: the
build's `bin/` and `scripts/` first on `PATH`, `SIMPLE_PATH` the build,
`SIMPLE_QSYS=local`, `SIMPLE_EMAIL`, `GFORTRAN_UNBUFFERED_ALL=y`).
`simple_test_exec` fixes its own seed; the `simple_exec` cases set
`SIMPLE_SEED`. The distributed master takes its thread count from OpenMP,
not from `nthr`, so `OMP_NUM_THREADS=1` is required.

| case | how | what it covers |
| --- | --- | --- |
| SW1 | `simple_test_exec test=simulated_workflow suite=1jxy` | the whole simulated workflow on about 100 particles: `solve2D` and its 2D table, single-state `prob` and `prob_neigh`, the nonuniform filter, shared memory; and its step 11, `solve3D` with `nparts=2`: the reprojection model, `prob_align` over two part files, volume assembly |
| SW6 | `... suite=6vxx` | as SW1 in shared memory, with the symmetry search |
| MSI | `simple_exec prg=refine3D projfile=<fixture>/7_solve3D/<project> nstates=3 refine=prob objfun=euclid nspace=1000 maxits=2 filt_mode=nonuniform nparts=2 nthr=1`, with the suite's point group (`d2`; 1jxy is beta-galactosidase) and mask diameter, no start volumes; then `refine3D refine=prob_neigh prob_neigh_mode=state nspace=2500 maxits=2` with the same keys on its output project and its three maps as `vol1` to `vol3` | the dense three-state table with its shift refinement and grid fit, its two part files, `ref_assign` with several states, the multi-state volume assembly and filter; the neighbourhood table in state mode |
| CCI | from the same project: `refine3D refine=prob objfun=cc nspace=1000 maxits=2 nparts=2 nthr=1`, plus `update_frac=0.5 trail_rec=yes` if `refine3D` takes them | the correlation kernels, which are separate from the Euclidean ones; fractional sampling and trailing reconstruction, if taken |
| PCI | as CCI with `objfun=euclid rec_backend=pcg filt_mode=nonuniform` | the PCG backend's call of the nonuniform filter |
| NU3D | `simple_exec prg=nu_filt3D` and `prg=postprocess_nu` on the fixture's final half maps | the filter's standalone programs, which no workflow case reaches |
| AD1 | `simple_test_exec test=solve3D_addon nthr=1` | 3,000 particles, fractional sampling, trailing reconstruction, `prob` and `prob_neigh`, the filter |
| CR1 | `simple_test_exec test=cont_refine3D_1jxy nthr=1` | `refine=cont`, `refine3D_auto` with `prob_neigh`, the filter and trailing reconstruction |

refine3D defaults to `mkdir=yes`, so the fixture's project is read, never
written. If a key of these command lines is refused, the case runs without
it and the journal says so.

- **MSI's start.** With `nparts` set and one state in the project, the
  distributed strategy draws uniform state labels and reconstructs the three
  start maps (`simple_refine3D_strategy.f90`, about lines 472-523; it stops
  if a state has five particles or fewer). SW1 has about 100 particles, so
  about 33 per state. If a state of MSI falls below 20 particles, build MSI
  instead on a larger simulated project (`simulate_particles` from the
  suite's truth map, about 1,500 particles, the suite's settings, true poses
  as the start), and record the change.
- **The fixture.** One base run of SW1 is made in a fixed directory that is
  never moved or rerun: MSI, CCI, PCI and NU3D read its outputs, and its
  project points at its stacks by absolute path.
- **The same paths.** Each case runs in one fixed directory for every build;
  after a run its outputs move to a directory named by build and case. With
  the same paths, files that embed absolute paths compare equal. AD1 and CR1
  work in a directory named by their process id, so for them only the check
  values in the log are compared.
- **What is compared.** Every regular file in a case's output tree, with
  exclusions per file or per field, all declared in Phase 0 in a list kept in
  the run directory, each with its reason. An exclusion is never added
  because an item changed a file. The known ones:
  - project files (`*.simple`) record the date and time of every job
    (`append_job_descr2jobproc`), so they are compared through text dumps
    instead: `simple_exec prg=print_project_field` of `ptcl2D`, `ptcl3D`,
    `cls2D`, `cls3D` and `out`, written outside the case tree; the dumps are
    compared outputs;
  - `solve3D_manifest.txt` without its run identifier and digest lines;
  - `sigma2_state.bin` without its layout digest (bytes 121-128), if the
    fixed paths do not already make it equal;
  - logs and `BENCH/` reports, which hold timings; a log is compared only
    after removing declared line patterns, where it carries results.
  The run `test_truth` arrived at the same rules
  (`simple_test_identity.f90` in its diff); its tool is not on master, so
  this run writes its own comparison script under `~/speed_mac`.
- **Determinism first.** In Phase 0 each case runs twice on the base build.
  A case whose compared outputs differ between those two runs cannot judge
  identity; the cause is found and fixed before any item.
- **Parallel code.** At one thread an OpenMP loop runs serially, so the
  identity cases cannot show that parallel code is safe. Every loop this run
  parallelises (the reads of the part files) and every per-thread buffer it
  adds (the reference copy, the repeated objective, the score-loop buffers)
  gets a check in its tester that runs the routine at one and at four threads
  over many inputs and compares the results bit for bit.
- **When to run what.** After every item: the item's build, the fast-gate
  suites, SW1, MSI, CCI, PCI and NU3D. At the end of every phase, in
  addition: SW6, the full Debug gate, and AD1 and CR1. AD1 and CR1 take an
  hour or more at one thread; they run in the background while the next
  phase starts, and a failure points at the phase they ran after.

## The items

**Phase 1: the probability tables** (`simple_eul_prob_tab.f90`,
`simple_eul_prob_tab_neigh.f90`, `simple_eul_prob_tab_utils.f90`,
`simple_commanders_prob.f90`, `simple_eul_prob_tab_tester.f90`)

1. The neighbourhood table read in one pass, with the counts written by the
   part (above), or dropped.
2. The compact record for the dense table: its own entry type without
   default values, writer, reader and empty-part writer, a format word in the
   header; the master holds the rows and the refined lists; the tester's
   `set_dist` follows; the dense table's double initialisation goes.
3. The parallel read of the dense table's part files, with the coverage check
   (every particle in exactly one part), the seed-table check after the loop,
   and the tester's comparison at one and several threads.

**Phase 2: the polar kernels** (`src/main/pftc/`)

4. The shift-matrix conversion limited to the searched shells.
5. The repeated objective after L-BFGS-B reused.
6. One double-precision copy of the reference per minimisation; the dead
   double-copy path and `shift_ref` deleted.
7. The score loop's per-reference allocations and the internal distance
   function replaced by per-thread buffers and a module routine, if the
   compiler confirms a trampoline or the allocations.

For items 4 to 7 the in-plane tester gets exact checks (`==`, not a
tolerance): an evaluation with the loaded reference against one that copies;
the trimmed conversion against the full one on the searched shells; the
searches over many reference and particle pairs at one and at four threads.

**Phase 3: volume assembly and the nonuniform filter**
(`simple_commanders_rec_distr.f90`, `src/main/nu_filt/`,
`simple_nu_state_filter.f90`, `simple_commanders_postprocess_nu.f90`,
`simple_commanders_resolest.f90`)

8. Timers: inside the filter (filtering, objective, optimiser, composition,
   writes) and per state, as lines of the volume-assembly report outside its
   "accounted for" total; the summary script `scripts/plot_refine3d_bench.py`
   shows them.
9. The two dead volume reads per state deleted; the filter's source halves
   taken from memory.
10. The filter's objective right after each candidate, without the cache
    write and read.
11. The composition filtering again from the input pair, kept only if the
    timing below says it is not slower; `postprocess_nu` and `resolest`
    follow; the cache files and their clean-up go once nothing writes them.

**Phase 4: the master's project handling** (`simple_refine3D_strategy.f90`,
`simple_commanders_prob.f90`, `simple_commanders_rec_distr.f90`,
`simple_commanders_euclid.f90`)

12. Volume assembly's trailing reread dropped; `refresh_matching_lp_from_project`
    reads one record.
13. The persistent reprojection builder, with the master reading the project
    first and the builder's project a copy of it (or a read of only the
    segments it needs).
14. The in-memory project handed to `prob_align`, volume assembly and
    `commit_sigma2` as copies, if the timing below shows that a copy is
    clearly faster than a read; otherwise items 12 and 13 stand alone and
    the measurement is reported.

**Timings that decide an item.** Two items are kept only on a measurement,
made with nothing else running on the Mac, three times each way:

- Item 11: `nu_filt3D` on one even/odd pair at box 256 or more, before and
  after. A pair can come from any half maps found read-only under
  `~/Processing`, or be made from the suite's truth map at a larger box with
  noise added independently to each half.
- Item 14: reading the largest project found read-only under `~/Processing`
  against copying it in memory (`sp_project` copy), in a throwaway program
  under `~/speed_mac/tools` compiled against the Release build's modules and
  library. No repository change for the timing.

## Common rules

- **Builds.** Release builds out of tree for the identity cases, each in its
  own directory with its own install (the install prefix defaults to the
  build directory), configured as `~/src/SIMPLE/build/CMakeCache.txt`
  (gfortran 16 from Homebrew, Release, `BUILD_TESTS=ON`, OpenMP,
  `USE_ARCHOPT=ON`). The base build is kept for the whole run. The Debug
  build lives in the worktree's `build/`: `./compile_debug.sh` (a clean
  build with the fast gate, whose 60-second budget only holds on an idle
  machine) runs once in Phase 0 and at each phase's end with nothing else
  running; per item the Debug build is incremental (`cmake --build build`)
  and the fast-gate suites run with `ctest -L "fast|provisional"` in
  `build/`, passing when no entry fails.
- **Gates per item:** the fast-gate suites pass on the incremental Debug
  build; `python3 scripts/check_test_registry.py .` passes whenever a tester
  changes; the identity cases are identical to the base.
- **No new tests or registrations.** New checks go into the existing tester
  modules (the 3D table tester, the in-plane tester, the filter's tester), so
  the test registry and the files the run `test_truth` edits stay untouched.
  Those files are: `production/CMakeLists.txt`, `production/simple_test_exec.f90`,
  everything under `src/main/commanders/test/` and `src/main/ui/simple_test/`,
  `simple_parameters.f90`, `simple_parameters_parse.f90`,
  `simple_commanders_sim.f90`, `simple_ui_sim.f90`, `simple_test_utils.f90`,
  `simple_forked_process.f90`, `simple_forked_process_tester.f90`,
  `simple_parameters_phases.f90`, `src/main/exec/simple_test_exec_highlevel.f90`,
  `src/main/apis/simple_test_exec_api.f90`, `scripts/simple_args_generator.pl`,
  `scripts/check_test_registry.py`,
  `src/main/flex/run/simple_flex_pca_application_tester.f90`,
  `doc/policies/test_environment_policy.md`,
  `doc/user_guide/simulate_particles.md` and
  `test_suite_truth_and_robustness.md`. A fix that needs one of them is
  recorded instead, for the maintainer. `test_truth` may also repair
  production code this run edits (the state handling of `ref_assign`, an
  empty state in the distributed PCG): `ref_assign`'s algorithm stays as it
  is apart from how it reads the distances, and the report lists every file
  both runs touched.
- **Fix as you go.** A defect found on the way is fixed in the phase that
  finds it, as an item of its own with its own diff. If the fix changes
  results, SW1 must pass its truth checks, CR1 its checks, and MSI must run
  to the end with every state populated; SW6 is reported but does not
  decide (its 80 picked particles make it pass on only some seeds). Then the
  base outputs are made again from a build with the fix, so the following
  items are still judged by identity. A design question is written up for
  the maintainer and the run goes on.
- **An item that cannot be made identical.** Find the cause. If it is the
  item's own mistake, fix it. If the item cannot be identical (for example
  the compiler vectorises a rewritten loop differently), restore its files
  from the previous snapshot with a plain copy (no hard links, no preserved
  times) and touch them, so that make rebuilds them; record the cause with
  the evidence, and go on with the next item. Never widen the comparison,
  never add an exclusion, never rerun until a case passes.
- **One diff per item.** After each item the source directories of the
  worktree (`src`, `doc`, `scripts`, `production`) are snapshotted with
  hard links (`rsync -a --link-dest`), and the item's diff is the
  difference between consecutive snapshots; the cumulative diff against
  HEAD is written too.
- **Nothing is committed**; the maintainer reviews the diffs.
- **CPU, memory, disk.** At most nine single-thread jobs at a time, the
  long ones started first (the Mac's cores are partly efficiency cores);
  builds take the cores left; the timings of items 11 and 14 run alone. The
  jobs' peak memories (measured in Phase 0 with `/usr/bin/time -l`) add up to
  at most 70% of the Mac's memory. An after-run whose outputs match is
  deleted once its dumps and logs are kept; the base outputs stay.
- **Not in this run.** Results-changing items; the worker fold; the parallel
  state assembly; the generated indexes under `doc/code_overview/`; skill
  files (proposed text goes into the report).

## Phases

**Phase 0: base and fixtures.** Build the worktree as it is (`./compile_debug.sh`
with the fast gate, on the idle Mac, and the base Release build). Run SW1
three times (once in the fixture directory first, then twice in its case
directory), then SW6, MSI, CCI, PCI and NU3D twice each, and AD1 and CR1
twice each, all on the base build, recording each job's wall time and peak
memory. Write the comparison script and the exclusion list from the pairs.
Exit: the gate passes; every case is deterministic under the list; MSI
keeps at least 20 particles per state (or the larger project is in place).

**Phases 1 to 4** carry out their items in order, with the per-item gates
and, at each phase's end, SW6, the full Debug gate, AD1 and CR1. Exit of
each: every kept item identical on every case, its diff written; items not
kept recorded with the cause. AD1 and CR1 of a phase may finish during the
next phase; a phase whose AD1 or CR1 differs is reopened before the run
goes on.

**Phase 5: documents and close.**

- `doc/policies/3D/refine3D_policy.md`: the dense table's record, the reads,
  the new benchmark lines; `doc/policies/NU/nonuniform_filtering_policy.md`:
  the filter without its disk cache, if item 11 is kept; the policy that
  describes volume assembly's report, if one does.
- `single_state_speed_week_completed.md`: the status of every item this run carried
  out, and the parallel state assembly reclassified as not identical, with
  the reasons above. `multistate_speed_run_completed.md`: one line that its item 3 was
  carried out here.
- The progress section of this note, then the report beside it,
  `identical_speed_run_report.md`: per item what changed, the identity
  evidence (cases, logs), the one-thread benchmark lines before and after as
  an indication only, the timings of items 11 and 14, what was not kept and
  why, the timing cases for the test machine, and a commit message (subject,
  blank line, short body).

**Phase 6: the rulings of 2026-10-18** (added after the report). Volume
assembly hands its project on to `commit_sigma2` and to the master, and
takes the master's project where it equals the file (ruling 3, identical);
then the label loss in volume assembly is fixed (ruling 2, a fix that
changes results on the path that back-fills projection labels). Each is
its own item. The only gate, by ruling: the Debug build and the fast-gate
unit suites; the Dell tests the items.

**Timing later, on the test machine.** After `test_truth`. First SW1, MSI
and CCI at one thread on the base and the final build there, since identity
shown with this Mac's ARM build need not carry over to the x86 build. Then
the cases of the runs `prob_speed` and `ms_speed` before and after this
run's diff, twice each way, interleaved: A and A-shm (`solve3D`,
beta-galactosidase), MS (three states from the class averages), M (the dense
three-state table at 1,000 directions), and E (`solve2D`).

## File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/completed/identical_speed_run.md` | this plan and its progress | each |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the compact record, the parallel read, the initialisation | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the one-pass read | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the dense entry type, `append_or_replace`; the sampler's buffer | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | `set_dist`, the read at one and several threads | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the reads, the empty-part writer; the project hand-off | 1 4 |
| `src/main/pftc/simple_polarft_geom.f90` | the conversion range; `shift_ref` deleted | 2 |
| `src/main/pftc/simple_polarft_corr.f90`, `simple_polarft_calc.f90`, `simple_polarft_core.f90` | the reference buffer, the dead path, the score-loop buffers | 2 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | the load at `minimize`, the repeated objective | 2 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | the exact checks | 2 |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90` | timers, dead reads, source halves; the project hand-off and the trailing reread | 3 4 |
| `src/main/nu_filt/simple_nu_filter_bank.f90`, `simple_nu_filter_state.f90`, `simple_nu_filter_apply.f90`, `simple_nu_filter.f90` | the filter without its cache; timers | 3 |
| `src/main/nu_filt/simple_nu_filter_tester.f90` | checks of the filter without its cache | 3 |
| `src/main/volume/simple_nu_state_filter.f90` | the kills after the composition | 3 |
| `src/main/commanders/simple/simple_commanders_postprocess_nu.f90`, `simple_commanders_resolest.f90` | the cache's other users | 3 |
| `scripts/plot_refine3d_bench.py` | the new benchmark lines | 3 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the persistent builder, the read order, the one-record read | 4 |
| `src/main/commanders/simple/simple_commanders_euclid.f90` | `commit_sigma2` on a copy | 4 |
| `doc/policies/3D/refine3D_policy.md`, `doc/policies/NU/nonuniform_filtering_policy.md` | the living documents | 5 |
| `doc/implementation_notes/completed/single_state_speed_week_completed.md`, `multistate_speed_run_completed.md` | status | 5 |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90`, `simple_commanders_euclid.f90` | the project handed on; the label fix | 6 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the hand-offs in the distributed master | 6 |

## Progress

### Phase 0, 2026-10-17: base and fixtures

No source file changed. The base is master `25d256886`; scratch is
`~/speed_mac`, the run directory `~/agent_runs/runs/speed_mac`.

- **Builds.** `./compile_debug.sh` passed its fast gate, 15 of 15 suites
  in 12.9 s. The base Release build is `~/speed_mac/builds/base`
  (gfortran 16.2 from Homebrew, Release, `BUILD_TESTS=ON`, OpenMP,
  `USE_ARCHOPT=ON`, installed into the build directory).
- **Cases.** `~/speed_mac/bin/run_case.sh` runs a case in its fixed
  directory with the CTest environment, `OMP_NUM_THREADS=1` and
  `SIMPLE_SEED=20261017`, dumps every project with `print_project_field`
  and moves the outputs aside. CCI and PCI took `update_frac=0.5
  trail_rec=yes`. NU3D filters copies of the half maps inside its own
  tree, since `nu_filt3D` writes beside its inputs (its first pair wrote
  into the fixture; those two files were deleted and the fixture checked
  against its checksums).
- **Determinism.** Every case ran twice on the base build with no
  differences under the declared exclusions (`exclusions.tsv` in the run
  directory): SW1 689 files, SW6 483, MSI 226, CCI 50, PCI 53, NU3D 14,
  AD1 44 check values, CR1 59.
- **Exclusions.** Project files (compared through their dumps), `BENCH/`,
  the time and process-id fields of `metadata.json`, the qsys job records
  (`EXIT_CODE_JOB_*.job`), the particle previews named by process id, the
  run identifier and checksum of `solve3D_manifest.txt`, and the timing
  lines of the PCG diagnostics. `sigma2_state.bin` needed none.
- **MSI** ends with 20, 30 and 37 particles per state, so the fixture's 87
  particles serve; the larger project was not needed.
- **Correctness checks on this Mac, one thread.** SW1 passes both steps.
  SW6 fails its volume-correlation floor (0.780 against 0.80), which the
  plan expects. AD1 fails `cohort_pose_error_excess_deg` (6.00 against
  5.0) in both runs; it is the floor the run `ms_speed` found sensitive to
  the draw. CR1 passes.
- **Cost** (wall time, peak memory): SW1 16 min, 3.2 GB; SW6 10 min,
  3.6 GB; MSI 5.5 min, 3.6 GB; CCI 0.5 min, 1.6 GB; PCI 3 min, 4.2 GB;
  NU3D 1 min, 0.6 GB; AD1 40 min, 1.7 GB; CR1 53 min, 10.2 GB. CR1 runs
  with at most AD1 beside it.

### Phase 1, 2026-10-17: the probability tables

All three items are kept; each is byte-identical to the base on SW1, MSI,
CCI, PCI and NU3D, and SW6 is identical after the phase.

- **The neighbourhood read in one pass.** The part writes the number of
  candidates of each particle after the particle indices. The master's
  first pass reads only these counts (and checks their sum against the
  header) instead of scanning every chunk; the second pass is unchanged,
  so candidates keep their file order. Simpler than before: kept.
- **The compact record for the dense table.**
  - `prob_dense_entry` (distance, in-plane index; no default values) and
    `prob_refined_entry` (particle, reference, shift) in
    `simple_eul_prob_tab_utils.f90`; a format word (-2) in the header slot
    where a sparse file holds its candidate count.
  - The worker keeps a batch's rows and refined entries and writes them at
    the particles' positions; `write_tab` checks that every particle's row
    was written. The master holds the rows and the refined entries by
    particle; `dense_candidate` rebuilds the entry `assign_candidate`
    reads. The table is no longer initialised before the read.
  - The sparse writer, its buffers and counters moved to the
    neighbourhood type; `write_empty_prob_tab` takes the caller's table,
    so each empty part keeps its own format. The unused `fill_tab` of both
    types is deleted.
  - On CCI's part file: 532 KB before, 153 KB after (19 particles, 1,000
    references, 17 refined entries).
- **The parallel read.** The dense part files are read in an OpenMP loop.
  After it, every particle of the master must be in exactly one part, and
  the parts must agree on the rotation grid (both checks are new).
- **Tests** (3D table tester, unit suite `unit_pftc_align2D3D`): a round
  trip of two part files written in batches and read by a master in
  reverse particle order, every field exact; and the read at one and at
  four threads over 1, 3 and 7 parts, bit for bit.
- **How the part files were compared.** They remain in the case trees and
  the items change their layout by design, so they are compared through a
  lossless decoder (every field, raw bits, per particle) instead of byte
  for byte; every other output is compared as declared in Phase 0.
- One-thread read times on these cases are hundredths of a second before
  and after; the gain shows only at scale (test machine).

### Phase 2, 2026-10-17: the polar kernels

All four items are kept; each is byte-identical to the base on SW1, MSI,
CCI, PCI and NU3D (each item built from its own snapshot), and SW6 is
identical after the phase.

- **The shift-matrix conversion.** `gen_shmat4aln` converts only the
  searched shells to single precision; all its consumers loop over
  `kfromto`.
- **The repeated objective.** The shift search keeps the last
  all-rotation objective it computed inside `minimize`, with its shift,
  and reuses it when the shift is bitwise equal: the L-BFGS-B callback
  scores every accepted point, and the alternation loop then scores the
  point the minimiser ended on. Outside `minimize` nothing is cached, since
  the references may change between calls.
- **One reference copy per minimisation.** A per-thread double-precision
  buffer (`pft_ref_aln_8`) is loaded at the start of `minimize`;
  `gen_corr_for_rot_8_aln` and `gen_corr_grad_for_rot_8_aln` call the
  unchanged kernels with it. The coarse searches, now private, use it too.
  `gen_corr_grad_only_for_rot_8` (unreachable) and `shift_ref` (no caller)
  are deleted.
- **The score loop.** The compiler confirmed both suspicions: heap
  trampolines for the internal distance functions
  (`__gcc_nested_func_ptr_created` per call; `-Wtrampolines` does not
  report heap trampolines) and a `malloc` of the score vector per call in
  `gen_likelihood_val_bank`. The likelihood routes now fill per-thread
  buffers and call `sample_likelihood_dists`, which shares its draw with
  `sample_likelihood_dist` (unchanged, still used by `ref_assign`). Two
  trampolines remain, in the bounded and power samplers used by the 2D
  sampling strategies and `greedy_sub`; they are outside the table fill
  and are left for a later change.
- **Tests** (in-plane tester and 3D table tester, all exact): the scores do
  not depend on the shells above the searched band; the loaded reference
  evaluates as the copied one, with other evaluations in between; the
  shift searches (24 references, three starts, both objectives) and the
  likelihood scores (shifted and bank routes) at one and at four threads;
  the array sampler draws as the function sampler on 200 random sets.
- **Indication only** (MSI, one thread, base run beside more jobs): the
  dense table's score loop 4.30 s before, 2.23 s after; its refinement
  0.30 s and 0.19 s.

### Phase 3, 2026-10-17: volume assembly and the nonuniform filter

All four items are kept; each is byte-identical to the base on SW1, MSI,
CCI, PCI and NU3D (each item built from its own snapshot), and SW6 is
identical after the phase.

- **Timers.** The filter times its candidate filtering, objective, label
  optimisation and composition (`get_nu_filter_timings`); the
  volume-assembly report gains a section "SPLIT TIMINGS" with these and,
  per state, the restore and the filter, outside the "accounted for"
  total. `scripts/plot_refine3d_bench.py` adds them as columns of
  `bench_iterations.csv`; its report lists them already. On MSI's
  three-state step the label optimisation is the largest part (8.7 s of
  14.9 s), the composition the smallest (0.2 s).
  - A first version freed the per-state timer arrays before the report
    read them; the three-state MSI case crashed on it (the single-state
    cases did not), and the item was corrected before its identity check.
- **The assembly's reads.** The two reads of the restored halves into
  `build%vol` and `build%vol2` are gone (nothing used them), and the
  filter's source halves are kept from the restoration instead of being
  read back from the files just written (the box only, as a read gives
  it).
- **The objective per candidate.** `setup_nu_dmats` filters each
  candidate pair from the transforms of the tapered inputs and scores it
  at once; the objective reads only the box, which the written file held
  exactly. The two transforms now live through the loop: two volumes more
  during the setup.
- **The composition without the cache.** `nu_filter_vols` always filters
  the source pair again per label (the competition pair, or the
  solvent-prior pair), with the same preparation and filter as the
  candidates; the base pair is released after the composition, and
  `nu_filt3D` passes its inputs. Nothing writes the cache any more, so its
  file names, writes, deletions and the unused `nu_filter_vol` are
  deleted. The filter tester checks that every composed voxel is its
  label's candidate, exactly.
  - The timing that decided it: `nu_filt3D` on a PfCRT pair (box 300,
    8 threads, alone on the Mac, interleaved): 10.56, 10.58 and 10.27 s
    with the cache, 9.95, 9.85 and 9.75 s without; the same peak memory
    (1.89 GB). On a slower or shared disk the gain should be larger.

### Phase 4, 2026-10-17: the master's project handling

All three items are kept, the third in part; each is byte-identical to the
base on SW1, MSI, CCI, PCI and NU3D, and SW6 is identical after the phase.

- **Fewer rereads.** Volume assembly no longer rereads the particle field
  for its trailing counts when the field in memory is the one it read at
  the start. When it has back-filled missing projection labels, it still
  rereads, as before; that reread drops the labels it has just filled in, an
  old inconsistency left for the maintainer. `refresh_matching_lp_from_project`
  reads the one particle record it uses through the binary reader:
  `sp_project%read_segment` ignores a record range for project files, so a
  range read through it reads every record.
- **The kept reprojection builder.** The distributed strategy keeps the
  builder of the reprojection model for the run. Its toolboxes are made on
  the first iteration (with the full read as before); later iterations read
  only the segments the model uses (the non-data segments, `out`, `optics`
  and the particle field), rebuild the small strategy toolbox and the three
  volumes the model's clean-up releases, and release the project after
  writing the field.
  - Two mistakes of a first version, both caught by the identity cases:
    the killed volumes were reused (a crash in the second iteration), and
    without `update_projinfo` the builder's project still named the project
    the run had been copied from, so the field was written into the
    fixture's project. The fixture was restored byte for byte (its Phase 0
    checksum), the runs that had read it were discarded, and every case now
    checks the fixture's checksums.
- **The project copy for `prob_align`.** Copying the largest project at hand
  (751,729 particles, 438 MB) takes 0.25-0.31 s, reading it 1.27-1.43 s.
  `prob_align` runs right after the master's full read, so the distributed
  master hands it a copy (`execute_from_project`; the builder takes an
  optional source project). Volume assembly and `commit_sigma2` are not
  changed: the file volume assembly reads follows writes the master's
  memory is not shown to match, and the one `commit_sigma2` reads holds
  volume assembly's resolution fields, which the master never takes. Handing
  them a project would need volume assembly to pass its own on.

### Phase 5, 2026-10-17: documents and close

- `doc/policies/3D/refine3D_policy.md`: the part files of the probability
  tables (the dense record, the parallel read and its checks, the
  neighbourhood counts), the shift search's reference copy and reused
  objective, the kept reprojection builder and the copy for `prob_align`,
  the in-memory source halves of the filter, the trailing counts, and the
  `SPLIT TIMINGS` section of the volume-assembly report.
- `doc/policies/NU/nonuniform_filtering_policy.md`: the filter without its
  scratch files (section 8).
- `single_state_speed_week_completed.md`: the status of every item, and the parallel
  state assembly marked as not identical with the reasons;
  `multistate_speed_run_completed.md`: its item 3 carried out here.
- The report: `identical_speed_run_report.md`.


### Phase 6, 2026-10-18: the rulings

By ruling, the only gate was the Debug build and the fast-gate unit suites
(15 of 15 pass after each item); no identity case was run. The Dell tests
both items.

- **Item 21, volume assembly hands its project on (ruling 3).** In the
  distributed gridding path (not `refine=sigma`), volume assembly takes a
  copy of the master's project (`execute_hand_on`) instead of reading the
  file: the master read it at the start of the iteration and
  `merge_algndocs` wrote the field from it. Assembly hands its project on
  as it ends; the low-pass refresh and `commit_sigma2`
  (`execute_from_project`) read that instead of the file. No write changes.
  The shared-memory path is unchanged. Meant to be identical.
- **Item 22, the back-filled labels are kept (ruling 2).** The trailing
  counts use the field in memory; the reread, which dropped the
  back-filled projection labels on trailing runs, is gone. A fix: results
  change on the path that back-fills.
