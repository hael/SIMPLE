# The 2D structure aligned with 3D, and one project owner on the master

Refactoring note, 2026-10-20. **Plan; nothing implemented.** The execution
plan of the unattended run `structure2D_owner` on the Oracle Linux test
machine (the Dell). The design comes from
`doc/implementation_notes/planned/speed_architecture_proposals.md`
(section 1.3, the row "one project owner on the master"; section 1.4 "The
commanders used as library routines"; section 1.5 "Moved before the
release", the 2D code structure) and from
`doc/implementation_notes/completed/identical_speed_run_review.md` ("The
master's project handling"). The code facts below were read at `34f314fe5`
on 2026-10-20; line numbers are approximate, find the code by name.

## Rulings

- 2026-10-19: the run takes the 2D structure first, then "one project owner
  on the master" for 3D and 2D, with no time limit. It starts by checking the
  last three speed changes on this machine's x86 build.
- 2026-10-20, random generator: each moved body builds its own `params`
  today, and `params%new` reseeds the random generator (with `SIMPLE_SEED`
  set, seed = base + 7919 x the number of earlier `seed_rnd` calls,
  `src/utils/math/simple_rnd.f90`) and sets the thread count. **No
  compatibility reseeds are kept**: once a body runs on the master's own
  `params`, the later random draws change. How those steps are judged is in
  "The tests".
- 2026-10-20, `cfar` and `lp_est`: when volume assembly works on the
  master's own project, the master sees `cfar` and `lp_est` of the current
  iteration instead of the previous one. Done as an isolated item after the
  rest, checked so that only `cfar`, `lp_est` and what reads them differ.
  The two items that need it (one write per iteration, and removing the last
  entry point) follow it.
- 2026-10-18 (the accepted proposals): the bodies of volume assembly and the
  other master steps move out of the commanders. The `simple-refine3d`
  skill's statement that volume assembly is intentionally executed in the
  commander is superseded; the report proposes the new skill text.
- Not in this run: the states' parallel assembly and the NU filter as an
  object (left alone, ruling of 2026-10-19), the worker fold, persistent
  workers, edits of stream code (`src/main/stream`, `commanders/stream`:
  Joe's area) and of the GUI (`nice/`).

## What the code is today

**3D.** `exec_refine3D` (`simple_commanders_refine3D.f90`) runs the iteration
loop and picks a strategy (`create_refine3D_strategy`,
`simple_refine3D_strategy.f90`): distributed when `nparts` is set and `part`
is not, shared memory otherwise. The matcher `refine3D_exec`
(`simple_strategy3D_matcher.f90`) samples, searches, writes the orientations
and sigma2 rows and writes the partial reconstructions of its particles. The
strategy assembles the volumes and decides convergence (`check_conv3D`) in
both modes.

**2D.** `exec_refine2D` (`simple_commanders_refine2D.f90`) runs the loop and
picks `create_refine2D_strategy` (`simple_refine2D_strategy.f90`) by the same
rule. The matcher `refine2D_exec` (`simple_strategy2D_matcher.f90`) samples,
searches and accumulates the class sums; then its
`finalize_restoration_and_convergence` either writes the part's contribution
(a distributed worker) or, in shared memory, restores the class averages
(`cavger_commit_carryover`, `cavger_restore_cavgs`, `cavger_gen2Dclassdoc`,
`write_merged`/`write_eo`), writes the `cls2D`/`cls3D` segments, runs
`check_conv2D`, applies `minits`, `maxits` and full-assignment coverage,
prints the verdict and writes the progress file. The distributed master
assembles through the `cavgassemble` commander
(`simple_commanders_mkcavgs.f90`, called from
`run_distributed_cavg_assembly`) and decides convergence in the strategy.
`refine2D_exec` is called by the shared-memory strategy and the worker
commander only. The stream's 2D pool always runs distributed
(`prg=refine2D_distr`, its workers `prg=refine2D`), never the shared-memory
branch.

**The master's project in distributed `refine3D`.** The strategy calls typed
entry points added to commanders: `commander_prob_align%execute_from_project`
(and the neighbourhood variant), `commander_volassemble%execute_hand_on` and
`commander_commit_sigma2%execute_from_project`. Their bodies
(`run_prob_align` and `run_prob_align_neigh` in `simple_commanders_prob.f90`,
about 90 lines each; `run_volassemble` in
`simple_commanders_rec_distr.f90`, about 530 lines, beside
`restore_state_from_parts` and `filter_pcg_nonuniform_maps`;
`run_commit_sigma2` in `simple_commanders_euclid.f90`, about 100 lines) each
build their own `params` and toolboxes from a command line and a copy of the
master's project. A kept reprojection builder (`refbuild`,
`materialize_reprojection_model_kept`) rereads the project and remakes
volumes by hand every iteration. Per iteration (gridding, `prob`, Euclidean
objective): two full project reads, four copies in memory, four writes of the
particle segment plus the `out` segment. The 2D master does the same through
`execute(cline)` of `prob_align2D`, `cavgassemble` (which writes the whole
project) and the sigma commit: three full reads, two segment rereads, one
full write and four segment writes per iteration.

**Other callers.** Volume assembly is also reached through
`simple_rec3D_service%assemble_gridding` (reconstruct3D, FLEX, the final
reconstruction, the solve3D add-on); the sigma commit through the solve3D
deferred commit and both 2D strategies; all of them through the private exec
programs. These keep working through the commanders' `execute(cline)`.

**Side effects of `params%new`** in each body (`simple_parameters_phases.f90`):
the reseed above, and `omp_set_num_threads(nthr)` with `nthr_glob`. A body
that runs on the master's `params` must set the thread count it needs itself
and the master must restore its own count after the call.

## The target

- **2D like 3D.** The 2D matcher searches, accumulates the class sums and
  writes orientations and sigma2 rows; a distributed worker also writes its
  contribution. Restoration of the class averages in shared memory and the
  convergence decision in both modes are strategy routines. `refine2D` is a
  public program, as `refine3D` is.
- **One project owner.** The bodies of `prob_align`, `prob_align_neigh`,
  volume assembly, the sigma commit, `prob_align2D` and class-average
  assembly become strategy (or domain) routines that take the master's
  `(params, build)` and explicit arguments for what differs today through the
  command line (threads, `which_iter`, `update_frac`, `combine_eo`, state),
  as `execute_rec3D_pcg_distributed_master(params, build, ...)` already does.
  Each commander's `execute(cline)` builds its `params` and builder and calls
  the same routine. The distributed master materialises its reprojection
  model on its own builder (the kept builder goes), reads only what the
  workers write, and writes the project once per iteration, at the end of the
  iteration, so the file always holds a complete iteration (the stream pool's
  hard stop and stage p07's `solve3D` rely on that). The typed entry points
  and the builder's `spproj_src` option are deleted.

## The items

### Phase 1: the 2D structure aligned with 3D

1. **Restoration and convergence out of the 2D matcher.** The matcher keeps
   everything through the orientation writes, the worker's
   `cavger_write_contribution` and `cavger_kill`, and
   `qsys_declare_part_finished`; in shared memory it returns with the sums in
   the class averager. New strategy routines: the shared-memory restoration
   (the matcher's present restoration lines, unchanged in order) and one
   convergence routine both modes call (`check_conv2D`, `minits`, `maxits`,
   coverage, the verdict, the progress file). The matcher's `converged`
   argument and the dead `refine2D_ctrl%display` go. Keep the order that
   makes results identical:
   - convergence runs before the shared-memory strategy restores
     `update_frac` (`check_conv2D` reads `params%l_update_frac`);
   - the carry-over flag (`l_frac`) is returned by the matcher (it is
     computed inside the matcher's set-up, and the shared-memory strategy
     may switch fractional update off before the matcher runs);
   - convergence runs after `cavger_gen2Dclassdoc`, which writes each
     particle's `res` that `check_conv2D` averages;
   - `check_conv2D` changes `params%trs` and `l_doshift`: shared memory
     applies them in the next iteration, distributed through the job
     description from iteration 4; keep both;
   - a fresh convergence object per iteration in shared memory, as today;
   - `cavger_new` stays in the matcher (reference preparation needs it);
   - shared memory is not routed through part files and `cavgassemble`: its
     particle field on disk is stale during the run.
2. **`refine2D` public.** A `new_refine2D` descriptor in
   `simple_ui_refine2D.f90` (category `refine2d`, the visibility `refine3D`
   has, required keys only `projfile` and `mskdiam`, since the private
   command of the same name is checked against the public required keys);
   `case('refine2D')` in `simple_exec_refine2D.f90` routed to
   `commander_refine2D`; the private `refine2D` stays the worker (the memory
   estimator keys on it); `scripts/default_audit.py` traces
   `set_refine2d_defaults`; a registration check in
   `simple_ui_visibility_tester.f90`; the UI JSON validator passes. The GUI
   accepts class selection only from `solve2D` (`nice/`): report it, do not
   edit.
3. **`solve2D_policy.md`** describes the new ownership (matcher, strategy,
   public program).

### Phase 2: the bodies out of the commanders (3D, pure moves)

4. **Volume assembly** (`run_volassemble`, `restore_state_from_parts`,
   `filter_pcg_nonuniform_maps`) moves into a new strategy module,
   `simple_volassemble_strategy.f90`, beside `simple_rec3D_pcg_strategy.f90`.
   The commander keeps building `params` and the builder and calls it; the
   strategy modules stop importing the commander module. Identical.
5. **`prob_align`, `prob_align_neigh` and the sigma commit** move the same
   way: the two alignment bodies into a new `simple_prob_align_strategy.f90`,
   the commit into `simple_strategy_sigma2.f90`. Identical.
6. **Volume names as outputs.** Assembly returns the names of the volumes it
   wrote instead of writing `params%vols`, `vols_even`, `vols_odd`, so that a
   caller's `params` never points at a state left without a map. Identical.

### Phase 3: the 3D master owns its project

7. **The master's own builder replaces the kept builder.** Build the general
   and 3D search toolboxes once on `build`; each iteration remake the volumes
   (a builder method, not hand-written in the strategy), rebuild the search
   toolbox, set the band and materialise the model on `build`, as the
   shared-memory path does. The full read after the model becomes a read of
   the particle segment while assembly still writes the file. Check `projinfo`
   and `out` in the dumps.
8. **`prob_align` and `prob_align_neigh` on the master's builder.** Sampling
   happens in place; the body sets the thread count it needs (the workers'
   `nthr`, as its command line said) and the master restores its own after.
   No reseed (ruling): later draws change.
9. **The sigma commit on the master's builder**, with explicit `which_iter`
   and thread count, the same way.
10. **The two `execute_from_project` entry points go.** The distributed
    strategy calls the bodies directly since items 8 and 9; shared memory
    stays on the commanders' `execute(cline)`. Identical to item 9's build.

### Phase 4: the 2D master owns its project

11. **`prob_align2D`, class-average assembly and the sigma commit for 2D**,
    in two steps:
    - *Pure move* (identity): the bodies of `exec_prob_align2D` and
      `exec_cavgassemble` become strategy routines (in
      `simple_prob_align_strategy.f90` and `simple_cavgassemble_strategy.f90`)
      that the commanders' `execute(cline)` call; the 2D strategy still calls
      the commanders. Identical to the Phase 3 build.
    - *On the master's project* (check build): the 2D strategy calls the
      routines with its own `(params, build)`. Class-average assembly writes
      only what it changes (`cls2D`, `cls3D`, `out` and the particle field) at
      the end of the iteration, not the whole project. The master's own
      `cls2D`, `out` and `jobproc` stay current in memory; its reread of the
      particle segment goes. Keep the 2D `update_frac` of the distributed
      master's `params` unchanged and pass `which_iter` as the commit needs it
      (`sigma2_commit_iter`). The `make_cavgs` callers keep working through the
      commander.

### Phase 5: assembly on the master's project, one write per iteration

12. **Volume assembly on the master's builder** (ruling: isolated). The low-pass
    refresh reads memory; `cfar` and `lp_est` reach the convergence log, its
    `LP_ESTIMATED` statistic and the final project one iteration earlier.
    Judged against the Phase 4 build: only `cfar`, `lp_est` and what reads them
    may differ. Any other difference: restore the item, record why, go on.
13. **One write per iteration.** A variant of `merge_algndocs` that merges into
    the master's field without writing; assembly and the PCG master stop
    writing the project; the segment rereads go. The master writes the
    particle segment once before the workers start (after the model and the
    sampling) and the iteration's results once at its end, plus `out`; the
    final full write as now. The distributed polish pass (`pose_cont=yes`)
    reproduces its sample from the project the merge wrote, so the merged
    particle segment is written before the polish pass starts. Identical to
    item 12's build, unless the Nyquist clamp of `res` in assembly ever binds
    (check: the clamp against `refresh_resolution_fields_from_fsc`).
14. **`execute_hand_on` and the builder's `spproj_src` option go**; the
    commanders' `execute(cline)` are thin wrappers.

### Phase 6: final checks, documents and close

15. **The truth suite on the final build.** `simple_test_exec test=nightly`
    on the base build and on the final build. Every entry that passes on the
    base passes on the final build; entries that already fail on the base (for
    example `solve3D_pyramid_blobs`, 19 of 20 in the last profile) are recorded
    as such. The workflow entries judge what identity cannot (the changed
    random draws of Phases 3 to 5).
16. **Documents.** `refine3D_policy.md` (the master's project, the bodies as
    strategy routines), `solve2D_policy.md`, the multi-state memory policy if
    the assembly's ownership text changes,
    `speed_architecture_proposals.md` (status), this note's progress, and the
    report beside the completed note, ending with a commit message.

## The tests

**Identity** (one thread, fixed seed `SIMPLE_SEED`, each case in one fixed
directory, outputs moved aside after each run, every regular file compared
except exclusions declared in Phase 0 with reasons, project files compared
through `print_project_field` dumps of ptcl2D, ptcl3D, cls2D, cls3D and out,
table part files through lossless decoders of their layouts). The
determinism of every case is shown in Phase 0 by two runs on the base build.
The run `test_truth_sync` left identity scripts and decoders on this machine
in `~/agent_runs/test_truth_sync/scratch/logs/identity/` (read, copy, never
write); the Mac runs' decoders know the newer 2D layouts (format words -4 and
-5) and are not on this machine: teach the copied decoder those layouts.

**Cases** (all made on the base build; the 3D cases read one fixture, the
SW1 directory of the base build, in every build):

| case | what it covers |
| --- | --- |
| SW1, SW6 | `simulated_workflow_1jxy` (with its two-part step) and `_6vxx`, including shared-memory `solve2D` |
| MSI | three-state distributed `refine3D` from the fixture, dense then neighbourhood table, `box_crop=128` so that the reconstruction's cropped records run |
| CCI | correlation objective, fractional update, trailing reconstruction (`box_crop=128`) |
| PCI | the PCG backend with the nonuniform filter (`box_crop=128`) |
| NU3D | `nu_filt3D` and `postprocess_nu` on the fixture's final half maps |
| SIG | distributed `refine3D refine=sigma` from the fixture |
| REC | distributed `reconstruct3D` from the fixture |
| MKC | distributed `make_cavgs` |
| E2 | `solve2D` with two parts on about 1,500 particles simulated from the fixture's truth map, all six stages (`snhc_smpl`, sparse `prob_snhc`, dense `prob`), cropping |
| E2S | E2's command line in shared memory |
| R2, R2S | `refine2D refine=prob` for two iterations from E2's final project, two parts and in shared memory |
| R2F | R2S with `update_frac` below 1 and carry-over, and `fillin=yes` |
| R2FD | R2F with two parts (distributed carry-over) |
| R2P | R2 through the public `simple_exec prg=refine2D` (from Phase 1 on; `mkdir=yes` puts the outputs in a subdirectory); its exclusions (the files that record the command line or the program name) are declared in Phase 0 by reason |
| CRD | distributed `refine3D` with the polish pass (`pose_cont=yes`) from the fixture, two parts: the polish pass reproduces its sample from the project the merge wrote |
| FLX | the weighted (state weight set) branch of volume assembly: the existing high-level test that reaches it (find it among the FLEX tests of `simple_test_exec`), compared by its check values |
| AD1, CR1 | `solve3D_addon` and `cont_refine3D_1jxy` at one thread, compared by their check values |

MKC runs distributed `make_cavgs` (two parts) on E2's final project. The
command lines of SW1, MSI, CCI, PCI and NU3D are in
`doc/implementation_notes/completed/identical_speed_run.md` ("The identity
test"), with `box_crop=128` added to MSI, CCI and PCI (MSI's second step takes
the first step's `recvol_state0N.mrc` maps); E2, E2S and R2 are described in
`doc/implementation_notes/completed/speed2d3d_run.md`. Each 2D case should
take at most about ten minutes at one thread.

**Judging each item.**
- Items that keep the random draws (Phases 1 and 2, items 6, 7, 10, 13, 14,
  and the pure-move step of item 11): byte-identical to the previous item's
  build on every case the item can reach.
- Items whose code runs a body on the master's `params` (8, 9, the second
  step of 11, and 12), so that the draws change by ruling: judge them on a
  **check build**: a copy of the working tree with the item applied plus a
  scratch-only patch that calls `seed_rnd` exactly where each `params%new` the
  item removed called it (the same calls, in the same order, together with the
  earlier seed items' patches). The check build must be byte-identical to the
  previous item's check build (or, for the first seed item, to its delivered
  build) on every case the item reaches; for item 12, apart from `cfar`,
  `lp_est` and what reads them (both are read only for reports). The patches
  never enter the repository; keep them in `scratch/keep/seed_patches/` and
  list them in the report. The delivered build then runs every case once to
  the end with finite results, and Phase 6 runs the truth suite.
- Items that build on others: 8 needs 7; 9 and 10 need 8; 12 needs 7 and 9;
  13 and 14 need 12. If an item is not kept, skip the items that need it and
  record why; the rest goes on.
- After every item: its Debug build and the fast gate, the unit suites of the
  area (`unit_pftc_align2D3D`, `unit_reconstruction`, the class and stream
  unit suites such as "pool 2D", "pool 2D object" and "class-average
  carry-over"), the identity cases it reaches; `scripts/check_test_registry.py`
  when tests or the test interface change. At the end of every phase: SW6,
  the full `./compile_debug.sh` on an idle machine, `mini_stream_*` and
  `lib_stream` CTest entries (the stream pool runs distributed `refine2D`),
  one multi-threaded E2S at 24 threads, MSI with two parts of 12 threads, and
  from Phase 4 on E2 with two parts of 12 threads (each to the end, finite).
  AD1 and CR1 at the end of Phases 1, 3 and 5.

**Parallel code.** No item adds a parallel region by design. If one does: no
string temporaries or module-buffer calls inside it, work arrays sized by the
team, errors raised after it, and a tester check at 24 threads repeated 20
times (the reason: commit 43de9b340).

## The x86 check of the speed changes (Phase 0)

The last three speed commits were shown byte-identical on the maintainer's
Mac (ARM) only. On this machine, at one thread, each against its parent:

| commit | parent | what must hold |
| --- | --- | --- |
| `aea5b1846` (finishing fixes, 2D speed) | `508d99066` | SW1, MSI, CCI, PCI, E2, E2S, R2 byte-identical; the 3D dense part files through the decoders of both layouts |
| `af13bc828` (3D reconstruction from the search's cropped records; cropped 2D class sums) | `aea5b1846` | MSI, CCI, PCI byte-identical (they run no 2D step; the 2D class sums change results by design) |
| `f078cf683` (2D tables, 3D table, NU label optimisation, 2D look and feel; the 2D table on the cropped image) | `0c9542e50` | MSI, CCI, PCI, NU3D byte-identical; E2, E2S and R2 with `cache=no` on `f078cf683` byte-identical to the same cases with `cache=yes` on `0c9542e50` (apart from the cache files and the files that record the command line); the 2D part files through decoders of both layouts |

Build each commit from an export (`git archive <commit> | tar -x -C <dir>`
under scratch); never move the repository's HEAD. Make each pair's inputs
(the fixture and E2's particle set and final project) with the pair's parent
build, so that no build reads a project written by a newer one. A difference is recorded
with the item of the speed run it belongs to and, if it is a defect, fixed in
Phase 0 as an item of its own; if it cannot be fixed identically, stop.

## Common rules

- The base is the HEAD the copy started from (`34f314fe5` or later) with the
  maintainer's uncommitted edits the copy carried.
- Release builds out of tree under `scratch/keep/build_*` with
  `BUILD_TESTS=ON`; the Debug build of `./compile_debug.sh` for the fast gate.
  `SIMPLE_PATH` and `PATH` set to the build in use in every run. After adding
  a source file, configure again (the source globs have no
  `CONFIGURE_DEPENDS`).
- Fortran as CLAUDE.md, the simple-modern-fortran skill and the compile-time
  policy say: `use ..., only:`, large types allocatable as components,
  `new`/`kill` symmetry, no `flag .and. allocatable_array(i)` on one line, no
  goto, THROW macros never continued with `//&`, no integer-overflow tricks,
  dead code deleted, short comments without history.
- Release 4 keeps no backwards compatibility except for reading `.simple`
  projects of earlier releases.
- Fix as you go: a defect found during the work is fixed in the phase that
  finds it, as an item of its own (add its file to the table first). A design
  question only the maintainer can answer is written up and the run goes on,
  unless it blocks the phase.
- Documents and the report are plain and self-contained: short sentences,
  terms defined, items named by what they do, the learning first.

## Phases

**Phase 0: base, cases and the x86 check.** Find how SIMPLE builds on this
machine; build the base (Release and Debug, fast gate) and the five commits of
the x86 check. Make the fixture (SW1 on the base build) and the E2 particle
set; run every case of "The tests" twice on the base build (AD1, CR1 and FLX
once); declare the exclusions from the pairs, and, by reason, the files the
R2P comparison and the x86 check's `cache=no` against `cache=yes`
comparisons may skip (files that record the command line or the program
name; the cache files). Run the x86 check. Exit: every case deterministic at
one thread under the exclusions; the 2D cases crop and reach their table
paths; MSI, CCI and PCI report the cropped records (`particle buffer +:
cropped records` in their refine3D reports); CRD runs its polish pass with
two parts; the x86 check's verdict recorded per commit.

**Phase 1** (items 1-3), **Phase 2** (items 4-6), **Phase 3** (items 7-10),
**Phase 4** (item 11), **Phase 5** (items 12-14): exit of each, every kept
item judged as "The tests" says, its evidence in the progress section; items
not kept recorded with the cause; the phase-end runs done.

**Phase 6** (items 15-16): the truth suite on the base and final builds, the
documents, the report, the note moved to `completed/`.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/refactoring_notes/planned/structure2D_and_project_owner.md` | this plan and its progress | each |
| `src/main/strategies/search/simple_strategy2D_matcher.f90` | items 1, 11 | 1 4 |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | items 1, 11 | 1 4 |
| `src/main/commanders/simple/simple_commanders_refine2D.f90` | items 1, 2 | 1 |
| `src/defs/simple_default_clines.f90` | item 2, if the defaults audit needs it | 1 |
| `src/main/simple_convergence.f90` | items 1, 12 | 1 5 |
| `src/main/ui/simple/simple_ui_refine2D.f90`, `src/main/ui/simple_ui_visibility_tester.f90`, `src/main/exec/simple_exec_refine2D.f90`, `scripts/default_audit.py` | item 2 | 1 |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90` | items 4, 6, 12, 14 | 2 5 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | items 5, 8, 10, 11 | 2 3 4 |
| `src/main/commanders/simple/simple_commanders_euclid.f90` | items 5, 9, 10, 11, 14 | 2 3 4 5 |
| `src/main/commanders/simple/simple_commanders_mkcavgs.f90` | item 11 | 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | items 4-10, 12-14 | 2 3 5 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90`, `simple_rec3D_pcg_strategy.f90`, `simple_make_cavgs_strategy.f90` | items 4, 6, 11, 13 | 2 4 5 |
| `src/main/strategies/parallelization/simple_strategy_sigma2.f90` | items 5, 9, 11 | 2 3 4 |
| `src/main/strategies/parallelization/simple_volassemble_strategy.f90` (new) | items 4, 6, 12-14 | 2 3 5 |
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` (new) | items 5, 8, 10, 11 | 2 3 4 |
| `src/main/strategies/parallelization/simple_cavgassemble_strategy.f90` (new) | item 11 | 4 |
| `src/main/strategies/search/simple_matcher_refvol_utils.f90` | item 7 | 3 |
| `src/main/simple_builder.f90` | items 7, 14 | 3 5 |
| `src/main/project/simple_sp_project_core.f90`, `src/main/project/simple_sp_project.f90` | items 11, 13 | 4 5 |
| `doc/policies/2D/solve2D_policy.md` | items 3, 16 | 1 6 |
| `doc/policies/3D/refine3D_policy.md`, `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`, `doc/policies/3D/reconstruct3D_pcg_policy.md` | item 16 | 6 |
| `doc/implementation_notes/planned/speed_architecture_proposals.md` | item 16 | 6 |

## Progress
