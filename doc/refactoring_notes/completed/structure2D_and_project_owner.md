# The 2D structure aligned with 3D, and one project owner on the master

Refactoring note, 2026-10-20. **Completed 2026-10-20** (all items kept; the
record is the "Progress" section, the summary the report
`structure2D_and_project_owner_report.md` beside this note). The execution
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
    base passes on the final build; entries that already fail on the base are
    recorded as such. The workflow entries judge what identity cannot (the changed
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
| MSI | three-state distributed `refine3D`, dense then neighbourhood table, `box_crop=128` so that the reconstruction's cropped records run; it starts from E2's 1,500 simulated particles with their simulated poses, because on x86 the fixture's 82 particles leave one state with a single particle after the first step and the second step refuses to start (Phase 0) |
| CCI | correlation objective, fractional update, trailing reconstruction (`box_crop=128`) |
| PCI | the PCG backend with the nonuniform filter (`box_crop=128`) |
| NU3D | `nu_filt3D` and `postprocess_nu` on the fixture's final half maps |
| SIG | distributed `refine3D refine=sigma` from the fixture |
| REC | distributed `reconstruct3D` from the fixture |
| MKC | distributed `make_cavgs` |
| E2 | `solve2D` with two parts on about 1,500 particles simulated from the fixture's truth map, all six stages (`snhc_smpl`, sparse `prob_snhc`, dense `prob`), cropping |
| E2S | E2's command line in shared memory |
| R2, R2S | `refine2D refine=prob` for two iterations from E2's final project, two parts and in shared memory |
| R2F | R2S with `update_frac=0.5` and `fillin=yes`, continuing from E2's last stage (`startit=2` with its carried class sums, `cavg_state.bin`), so that both iterations blend carry-over; a fresh start never blends carry-over (Phase 0) |
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
| `src/main/commanders/simple/simple_commanders_prob.f90` | the Phase 0 fix of the 2D sampler (its comment), items 5, 8, 10, 11 | 0 2 3 4 |
| `src/main/commanders/simple/simple_commanders_euclid.f90` | items 5, 9, 10, 11, 14 | 2 3 4 5 |
| `src/main/commanders/simple/simple_commanders_mkcavgs.f90` | item 11 | 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | items 4-10, 12-14 | 2 3 5 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90`, `simple_rec3D_pcg_strategy.f90`, `simple_make_cavgs_strategy.f90` | items 4, 6, 11, 13 | 2 4 5 |
| `src/main/strategies/parallelization/simple_strategy_sigma2.f90` | items 5, 9, 11 | 2 3 4 |
| `src/main/strategies/parallelization/simple_volassemble_strategy.f90` (new) | items 4, 6, 12-14 | 2 3 5 |
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` (new) | items 5, 8, 10, 11 | 2 3 4 |
| `src/main/strategies/parallelization/simple_cavgassemble_strategy.f90` (new) | item 11 | 4 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the Phase 0 fix: a fresh 2D start with `fillin=yes` samples by update count | 0 |
| `src/main/image/simple_accum_blend_tester.f90` | item 4 (its comments name the module that holds `blend_trailing_accumulators`) | 2 |
| `doc/policies/NU/nonuniform_filtering_policy.md`, `doc/policies/3D/solve3D_addon_policy.md` | item 16 (they name the volume-assembly commander as the owner of the moved code) | 6 |
| `src/main/strategies/search/simple_matcher_refvol_utils.f90` | item 7 | 3 |
| `src/main/simple_builder.f90` | items 7, 14 | 3 5 |
| `src/main/project/simple_sp_project_core.f90`, `src/main/project/simple_sp_project.f90` | items 11, 13 | 4 5 |
| `doc/policies/2D/solve2D_policy.md` | items 3, 16 | 1 6 |
| `doc/policies/3D/refine3D_policy.md`, `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`, `doc/policies/3D/reconstruct3D_pcg_policy.md` | item 16 | 6 |
| `doc/implementation_notes/planned/speed_architecture_proposals.md` | item 16 | 6 |
| `doc/refactoring_notes/completed/structure2D_and_project_owner.md`, `doc/refactoring_notes/completed/structure2D_and_project_owner_report.md` | item 16: this note moved to `completed/` at the close, and the report | 6 |

## Progress

### Phase 0, 2026-10-20: base, cases and the x86 check

**Verdict of the x86 check: all three speed commits are byte-identical to
their parents on this machine's x86 build.** One defect was found and fixed
(the 2D sampler below); two cases had to be defined differently from the
Mac's (MSI and R2F, in the table of "The tests").

Files changed: this plan;
`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` and
`src/main/commanders/simple/simple_commanders_prob.f90` (the sampler fix and
a comment that described the sampler wrongly). The run directory below is
`~/agent_runs/structure2D_owner` (`R`); scripts, decoders and logs are in
`R/scratch`, the builds, inputs and case trees in `R/scratch/keep`, the
journal in `R/state/journal_phase_0.md`.

- **How SIMPLE builds here.** GCC 15.2 (`gfortran` from `gcc-toolset-15`,
  first on `PATH`), CMake 3.26, OpenMP, `USE_ARCHOPT=ON`, no modules; the
  gate's Python scripts need Python 3.12 first on `PATH`
  (`R/scratch/pybin`). Release builds with `BUILD_TESTS=ON` out of tree
  (`R/scratch/build_release.sh`), about 1.5 minutes each at 24 jobs. The
  base is `abeeebfba` with no uncommitted edits. Builds kept:
  `keep/build_base` (the base) and `keep/build_p0` (the base with the
  sampler fix: the build later phases compare with). The Debug build of
  `./compile_debug.sh` passed the fast gate, 15 of 15, before and after the
  fix (`logs/compile_debug_base.log`, `logs/compile_debug_p0.log`; 107 s,
  over the 60 s budget only because cases were running).
- **The cases.** `R/scratch/run_case.sh` runs each case at one thread with
  `SIMPLE_SEED=20261019` in its fixed directory `keep/idcase/<case>/run`
  and moves the outputs to `out_<tag>`; projects are dumped with
  `print_project_field`; `R/scratch/compare_case.sh` compares two output
  trees. The inputs, all made once on the base build and never rerun: the
  fixture (`keep/idcase/fixture_base`, SW1); E2's particle set
  (`e2data_base`: 1,500 particles simulated from the fixture's 1JXY map,
  box 192, 1.3 A, SNR 0.1, shift error 3 pixels, defocus 2 +- 0.5 um,
  astigmatism 0.1, imported with only their CTF parameters); R2's start
  project (`r2data_base`: E2's command line run once); MSI's start project
  (`msidata_base`: the same particles imported with their simulated poses).
  The 3D cases with `box_crop=128`; R2P is defined but runs from Phase 1 on
  (`compare_r2p.sh`).
- **Determinism.** Every case ran twice on the base with no difference
  under the exclusions (files compared): SW1 641, SW6 581, MSI 226, CCI 50,
  PCI 53, NU3D 19, SIG 27, REC 29, CRD 51, MKC 19, E2 147, E2S 146, R2 32,
  R2S 30, R2F 30, R2FD 32. Run once: AD1 passes its 46 checks, CR1 its 60,
  FLX (`flex_pca_blobs`, the test that reaches the weighted state branch of
  volume assembly through its per-state `refine3D_auto` runs) its 95 rows
  of `metrics.tsv`. Logs: `logs/chain_<case>_base.log`.
- **What the cases reach.** E2 and E2S crop 192 to 88 and run stages 1-2
  with `snhc_smpl`, 3-5 with the sparse `prob_snhc` table and 6 with the
  dense `prob` table; R2 keeps the dense 2D part files (format word -5).
  MSI, CCI and PCI report `cropped records` in their refine3D reports. CRD
  runs the polish pass (continuous pose work) after each discrete pass, both
  with two parts. MSI's three states hold 448 to 554 particles.
- **The exclusions** (declared here, frozen from now on): project files are
  compared through their dumps; logs and `BENCH/` reports (timings); job
  records `EXIT_CODE_JOB_*`; particle previews named by process id;
  `metadata.json` without its time and process-id fields; PDF dates and
  document ids; the PCG diagnostics' timing lines; `solve3D_manifest.txt`
  without run id and checksum; table part files through lossless decoders
  (`R/scratch/decode/decode_tab.f90` reads all five layouts of the dense and
  2D part files, `decode_neigh.f90` the neighbourhood files; between two
  layouts only the fields both hold are compared). New on x86: a
  `metadata.json` is skipped where two processes published particle
  previews in its directory: in SW1's extraction step the master and its
  worker both write the graphical interface's metadata at the same moment
  and the last write wins (seen once in the x86 check; a race in the
  interface's metadata writer, `src/utils/gui`, reported for the
  maintainer). The `cache=no` against `cache=yes` comparison and the R2P
  comparison needed no file excluded: no output records the command line or
  the program name, and the cache files are gone at the end of a run
  (`R/scratch/exclusions.tsv`).
- **The x86 check** (each pair's inputs made on the parent build):
  `aea5b1846` against `508d99066`: SW1 642, MSI 226, CCI 50, PCI 53, E2
  147, E2S 146, R2 32 identical, the 3D dense part files compared across
  their two layouts (format word -2 against -3 with its active-state list).
  `af13bc828` against `aea5b1846`: MSI, CCI, PCI identical. `f078cf683`
  with `cache=no` against `0c9542e50` with `cache=yes` (in the same
  directory): MSI, CCI, PCI, NU3D, E2, E2S and R2 identical, R2's part files
  compared across the old candidate stream and the dense format -5; an
  extra one-iteration `prob_snhc` run with two parts (the sparse part files
  `solve2D` deletes) identical across the old stream and format -4, chunks
  and candidate counts included. Logs: `logs/chain_<case>_x1|x2|x3.log`.
- **Fixed: a fresh 2D start with `fillin=yes` never ended.** With
  `startit=1` the 2D sampler keeps the subset it drew first (the sticky
  stage of `solve2D`'s stage 1), while `fillin=yes` makes `refine2D` iterate
  until every particle has been assigned: half the particles were never
  sampled, and R2F as first written ran past 145 iterations. `solve2D` never
  combines the two. Now a fresh start with `fillin` samples by update count,
  as later stages do. Nothing that ended before changes: the fixed build is
  byte-identical to the base on SW1, SW6, E2, E2S, R2, R2S, R2F, R2FD and
  MKC (`logs/chain_<case>_p0.log`), and the fresh-start form (case R2F1)
  now ends after two iterations with every particle assigned.
- **MSI and R2F redefined.** MSI from the fixture fails on x86 (a state
  keeps one particle after the first step; the second step needs six), so
  it starts from E2's particles with their simulated poses, as
  `identical_speed_run.md` foresaw for a thin state. R2F as written never
  blended carry-over (blending needs a run that continues, `startit` above
  1) and hit the defect above; it now continues from E2's last stage with
  its carried class sums.
- **Cost** at one thread on a loaded machine (wall time, peak memory of the
  master): SW1 24 min, 1.9 GB; SW6 28 min, 2.0 GB; MSI about 10 min; PCI
  56 s; CRD 46 s; CCI, NU3D, SIG, REC under 15 s; E2 1.3 min; E2S 2.5 min;
  R2, R2F, R2FD under 20 s; R2S 33 s; AD1 25 min, 1.0 GB; CR1 34 min,
  8.5 GB; FLX 22 min, 7.2 GB.
- **References for Phase 1.** The case trees were moved under
  `keep/idcase` (the old `scratch/idcase` is a link to it, so that the
  inputs' absolute paths still resolve), and every case ran once more there
  on `build_p0` (tag `p0`): the trees Phase 1's first item is compared with.

Exit items: every case deterministic at one thread under the exclusions
(above); the 2D cases crop and reach their table paths; MSI, CCI and PCI
report the cropped records; CRD runs its polish pass with two parts; the
x86 check's verdict is recorded per commit (identical, all three).

### Phase 1, 2026-10-20: the 2D structure aligned with 3D

All three items kept; every case they reach is byte-identical to the build
before. Files changed: `simple_strategy2D_matcher.f90`,
`simple_refine2D_strategy.f90`, `simple_commanders_refine2D.f90` (item 1);
`simple_ui_refine2D.f90`, `simple_exec_refine2D.f90`,
`simple_ui_visibility_tester.f90`, `scripts/default_audit.py` (item 2);
`doc/policies/2D/solve2D_policy.md` (item 3); this plan. Builds
`R/scratch/keep/build_i1` and `build_i2`; logs in `R/scratch/logs`,
journal `R/state/journal_phase_1.md`.

- **Restoration and convergence out of the 2D matcher.** The matcher
  (`refine2D_exec`) now samples, searches, accumulates the class sums and
  writes the orientations and sigma2 rows; a distributed worker also writes
  its class-sum contribution; in shared memory it returns with the sums in
  the class averager and hands back whether the iteration blends carried
  sums. Its `converged` argument and the unused `refine2D_ctrl%display` are
  gone. Two new routines in the strategy: `restore_cavgs_inmem` (the
  matcher's former restoration lines, in their order) and
  `check_refine2D_convergence`, which both modes call (`check_conv2D`,
  `minits`, `maxits`, full-assignment coverage, the verdict, the progress
  file). Shared memory runs matcher, restoration, convergence (with a fresh
  convergence object per iteration, as before), then restores
  `update_frac`; the distributed master keeps its shift-search activation
  after the routine. The restoration now runs after the matcher has freed
  its polar transforms and sigma object and written its benchmark report;
  the restoration code reads neither. Identity against Phase 0's build:
  SW1 641 files, SW6 581, E2 147, E2S 146, R2 32, R2S 30, R2F 30, R2FD 32,
  R2F1 30, MKC 19, all identical (`logs/chain_<case>_i1.log`).
- **`refine2D` public.** `simple_exec prg=refine2D` runs in shared memory
  or, with `nparts`, as the distributed master; its descriptor (category
  "Refine 2D Workflows", the visibility `refine3D` has) requires only
  `projfile` and `mskdiam`, because the private worker of the same name is
  checked against the public required keys. The workers stay on
  `simple_private_exec` (the queue picks the executable, not the program
  name), and the memory estimator still keys on the worker's name. The
  defaults audit traces `set_refine2D_defaults`, so the descriptor shows
  `maxits` 30, `cls_init` `ptcl`, `refine` `snhc_smpl`. The UI JSON
  validator passes; `unit_ui` checks the category, the executable and the
  two required keys (894 assertions). Side effect: a command line with
  `prg=refine2D` run in-process now honours `mkdir` and project discovery;
  `solve2D` and `refine2D_nano` already set `mkdir=no` and the project.
  Identity against item 1's build: the same cases, all identical. R2P (the
  public route with two parts) equals R2 on 33 compared items once the
  execution directory is mapped in both of its forms: absolute, and
  relative to the project root as the star files write it (`1_refine2D/`
  where R2 writes `run/`); its case-root project is the input, which
  `mkdir=yes` leaves untouched. No file excluded.
- **`solve2D_policy.md`** describes the public program, the matcher's and
  the strategy's new shares, and the sampler change of Phase 0.
- **Checks after each item:** the Debug build with the fast gate, 15 of 15
  (`unit_pftc_align2D3D`, `unit_reconstruction`, and `unit_stream` with the
  "pool 2D", "pool 2D object" and "class-average carry-over" suites are in
  the gate).
- **Phase-end runs** on item 2's build: SW6 identical (above); `lib_stream`,
  `mini_stream_6vxx` and `mini_stream_1jxy` pass on the Debug build
  (`logs/ctest_stream_p1.log`); E2S at 24 threads and MSI with two parts of
  12 threads run to the end with no NaN or infinity (`logs/mt_i2.log`); AD1
  and CR1 pass with `metrics.tsv` identical to the base's; the full
  `./compile_debug.sh` on the idle machine passes, its gate 15 of 15 in
  16 s, within the 60 s budget (`logs/compile_debug_p1_final.log`).
Exit: every item kept and judged by identity, its evidence above; the
phase-end runs done.

- **For the report** (not edited here): the graphical interface accepts a
  class selection only from `solve2D` (`nice/nice_lite/views/job_builder_views.py`),
  so `refine2D` results cannot be selected there yet; `doc/user_guide/user_guide.md`
  does not list `refine2D`.

### Phase 2, 2026-10-20: the 3D bodies out of the commanders (pure moves)

All three items kept; every case they reach is byte-identical to the build
before. Files changed: the new
`src/main/strategies/parallelization/simple_volassemble_strategy.f90` and
`simple_prob_align_strategy.f90`; `simple_commanders_rec_distr.f90`,
`simple_commanders_prob.f90`, `simple_commanders_euclid.f90`,
`simple_strategy_sigma2.f90`, `simple_refine3D_strategy.f90`,
`simple_rec3D_service.f90`, `src/main/image/simple_accum_blend_tester.f90`
(comments only); this plan. Builds `R/scratch/keep/build_i4`, `build_i5`,
`build_i6`; journal `R/state/journal_phase_2.md`.

- **Volume assembly moved.** The routine that assembles the state volumes
  from the workers' partial reconstructions (`run_volassemble`, with
  `restore_state_from_parts`, `filter_pcg_nonuniform_maps` and their types
  and constants) is now in `simple_volassemble_strategy.f90`, unchanged. It
  still builds its own parameters and builder from its command line; the
  commander keeps its type with `execute` and `execute_hand_on` as one-line
  wrappers. Decision: the plan asks both that the commander keep building
  the parameters and that the strategy modules stop importing the
  commander; the refine3D strategy and the reconstruction service build
  nothing themselves, so the command-line entry moved with the body and
  both now call it from the strategy module. Running the body on a
  caller's own parameters and builder is the work of Phase 5. Identity
  against Phase 0's build (`logs/chain_<case>_i4.log`): SW1 641 files, SW6
  581, MSI 226, CCI 50, PCI 53, CRD 51, REC 29, SIG 27, NU3D 19.
- **Probability alignment and the sigma commit moved.** The dense and
  neighbourhood alignment bodies (`run_prob_align`,
  `run_prob_align_neigh`, with their output cleanup and benchmark writer)
  are in `simple_prob_align_strategy.f90`; the sigma commit
  (`run_commit_sigma2`) is in `simple_strategy_sigma2.f90`. Decision: in
  shared memory the alignment fills its table through the table commander,
  which lives in the module that calls the bodies, so importing it would be
  circular; the bodies take the table commander as an optional argument of
  the abstract commander type, which the alignment commanders pass and the
  distributed master, which always schedules parts, does not need. The
  strategies still call the commanders' entry points; Phase 3 replaces
  them. Identity against item 4's build: the 3D cases above and E2 147,
  E2S 146, R2 32, R2S 30, R2F 30, R2FD 32, R2F1 30, MKC 19 (the sigma
  commit is in every Euclidean run).
- **Volume names as outputs.** Assembly returns the names of the volumes
  each state ends with (restored, or carried forward from disk) in optional
  outputs instead of writing `params%vols`, `vols_even` and `vols_odd`, so
  that a caller's parameters never name a state left without a map; nothing
  read those fields after assembly wrote them. Identity against item 5's
  build: the 3D cases, E2 and R2.
- **FLX is not deterministic.** FLX (`flex_pca_blobs`, the weighted branch
  of volume assembly) runs its per-state refinements at two threads
  whatever `nthr` is, so its check values vary from run to run: two runs on
  item 4's build differ in 12 of its 95 rows, as the base run differs from
  item 2's build. Phase 0 ran it once and could not see this. It is judged
  by its check values passing their floors: no failure on the base, item
  2, item 4 (twice) or item 6.
- **Checks after each item:** the Debug build with the fast gate, 15 of 15
  for items 4, 5 and 6 (`logs/compile_debug_i4.log` to `_i6.log`, built
  in each item's snapshot).
- **Phase-end runs** on item 6's build: SW6 identical (above); E2S at 24
  threads and MSI with two parts of 12 threads ran to the end with no NaN
  or infinity, the master reporting 12 threads (`logs/mt_i6.log`); the full
  `./compile_debug.sh` on the idle machine passes, gate 15 of 15 in 15 s
  (`logs/compile_debug_p2_final.log`); `lib_stream` and the two
  `mini_stream` entries pass (`logs/ctest_stream_p2.log`).

Exit: every item kept and judged by identity, its evidence above; the
phase-end runs done.

### Phase 3, 2026-10-20: the 3D master owns its project

All four items kept. Files changed: `simple_refine3D_strategy.f90`,
`simple_builder.f90`, `simple_prob_align_strategy.f90`,
`simple_strategy_sigma2.f90`, `simple_commanders_prob.f90`,
`simple_commanders_euclid.f90`; this plan. Builds
`R/scratch/keep/build_i7` to `build_i10` (items as delivered) and
`build_c8`, `build_c9` (check builds); the seed patches and the thread
probe are in `R/scratch/keep/seed_patches`; journal
`R/state/journal_phase_3.md`.

- **The master's own builder replaces the kept builder.** The distributed
  refine3D master used to keep a second builder for its reprojection model
  (the references the workers search), reread the project into it each
  iteration, and then read the whole project file again into its own
  builder. Now the model is made on the master's builder: the general
  toolbox once, at the first iteration; on later ones the volumes again
  (`new_vols`, a new builder method that `build_general_tbox` also uses),
  the search toolbox, the band and the model. Before the model the master
  reads the segments the kept builder read (the non-data segments, `out`,
  `optics` and the particle field), because the particle field on file is
  volume assembly's version, and sets the project's file name as a full
  read does (`update_projinfo`); the full read after the model is gone.
  Identity against Phase 2's build: SW1 641 files, MSI 226, CCI 50, PCI 53,
  CRD 51, SIG 27, REC 29; the final projects' `projinfo` is unchanged.
- **An incident while testing it.** The first build of this item left out
  `update_projinfo`. A run started from a project in another directory
  then names that project in `projinfo`, so the master's writes went to
  the fixture's project (the shared input of the 3D cases). It was caught
  by the identity check (CCI differed). The jobs were stopped, and the
  fixture was regenerated on the base build at its original path. It
  equals the damaged copy in every file but that project, and the project
  equals SW1's reference run in every segment. The fixed build is the
  one above.
- **`prob_align` and `prob_align_neigh` on the master's builder.** The
  alignment bodies (`prob_align_exec`, `prob_align_neigh_exec`) take the
  caller's parameters, builder, command line and thread count; the
  commanders build their own and call them. The master draws the sample in
  its own particle field, runs the bodies at the workers' thread count and
  restores its own. The `params%new` of each body, which reseeded the random
  generator, is gone (ruling: no compatibility reseed). Check build (one
  `seed_rnd` call where that `params%new` ran, patch `p8.py`) against
  item 7's build: SW1, MSI, CCI, PCI, CRD, SIG identical. The delivered
  build ran every case it reaches to the end with finite results.
- **The sigma commit on the master's builder.** `commit_sigma2_exec` takes
  the caller's parameters and builder, an explicit `which_iter` and the
  thread count; the master calls it after volume assembly. Check build
  (patches `p8.py` and `p9.py`) against item 8's check build: SW1, MSI,
  CCI, PCI, CRD, SIG identical.
- **Thread counts.** A probe build (item 10 with log lines, `thread_probe.py`,
  never in the repository) ran MSI's first step with two parts of six
  threads while the master had twelve: in the second iteration the master
  stood at 12 before the alignment, the alignment ran at 6, and the master
  was back at 12 after it; the commit ran at the 12 its command line gives
  (`logs/thread_probe_t10.log`), finite.
- **The `execute_from_project` entry points gone** (of `prob_align`,
  `prob_align_neigh` and the sigma commit), with the project argument of
  their wrappers. Identity against item 9's build: SW1, MSI, CCI, PCI, CRD,
  SIG.
- **Checks after each item:** the Debug build with the fast gate, 15 of 15
  for items 7 to 10.
- **The 2D cases** reach the sigma commit's command-line wrapper, which
  items 9 and 10 changed: E2, E2S, R2, R2S, R2F, R2FD, R2F1 and MKC are
  identical on item 9's build against Phase 2's and on item 10's against
  item 9's.
- **Phase-end runs** on item 10's build: SW6 identical; AD1 passes with its
  check values identical to the base (it runs in shared memory); CR1
  passes its 60 checks, its values from stage e on differing from the base
  as its distributed refinements' random draws now do; E2S at 24 threads
  and MSI with two parts of 12 threads ran to the end with no NaN or
  infinity (`logs/mt_i10.log`); the full `./compile_debug.sh` on the idle
  machine passes, gate 15 of 15 in 15 s (`logs/compile_debug_p3_final.log`);
  `lib_stream`, `mini_stream_6vxx` and `mini_stream_1jxy` pass
  (`logs/ctest_stream_p3.log`).

Exit: every item kept and judged as "The tests" says (items 7 and 10 by
identity, items 8 and 9 on their check builds), its evidence above; the
phase-end runs done.

### Phase 4, 2026-10-20: the 2D master owns its project

Item 11 kept in both steps. Files changed: the new
`src/main/strategies/parallelization/simple_cavgassemble_strategy.f90`;
`simple_prob_align_strategy.f90`, `simple_refine2D_strategy.f90`,
`simple_commanders_prob.f90`, `simple_commanders_mkcavgs.f90`; this plan.
Builds `R/scratch/keep/build_i11a` (the pure move), `build_i11b` (as
delivered) and `build_c11` (its check build); seed patch `p11.py` and the
probe `thread_probe2D.py` in `R/scratch/keep/seed_patches`; journal
`R/state/journal_phase_4.md`.

- **The pure move.** The body of `prob_align2D` (sampling, the particle
  cache, the 2D probability table filled in process or on the parts, the
  class assignment) became `prob_align2D_exec` in
  `simple_prob_align_strategy.f90`, the body of class-average assembly
  (reducing the workers' partial class sums, restoring and writing the
  class averages, the class document) `cavgassemble_exec` in the new
  `simple_cavgassemble_strategy.f90`; each has a command-line wrapper that
  builds its own parameters and builder, which the commanders call. Identity
  against Phase 3's build: E2 147 files, E2S 146, R2 32, R2S 30, R2F 30, R2FD
  32, R2F1 30, MKC 19, SW1 657, SW6 581.
- **On the master's project.** The distributed refine2D master calls the
  three bodies (alignment, assembly, sigma commit) with its own parameters
  and builder: the alignment draws its sample in the master's particle field
  (a full update sets the master's `update_frac` to 1 for the call and
  restores it), the commit gets `which_iter` as `sigma2_commit_iter` gives
  it. Assembly updates the master's class document, `out` segment and the
  particles' resolutions in memory and the master writes only those
  segments (`cls2D`, `cls3D`, `out`, the particle field) instead of the
  whole project; its rereads of the particle segment after the alignment
  and after assembly are gone. The command-line route (`make_cavgs`, the
  private `cavgassemble`) keeps its full write, which also keeps the file's
  `projinfo` current. Shared-memory refine2D stays on the commanders, as 3D
  does.
- **A difference found and kept identical.** The first check build differed
  in the even and odd class averages of distributed runs: restoration picks
  the even/odd low-resolution averaging limit from the builder's previous
  class FRCs (Fourier ring correlations) before it stores the new ones, and
  a freshly built toolbox has none, while the master's builder kept the
  previous iteration's. Assembly now starts the class FRCs afresh, as a new
  toolbox did. For the maintainer: shared-memory refine2D reads the previous
  FRCs from file for the same averaging, so the two modes differ there.
- **Judged on the check build** (one `seed_rnd` call before each of the
  three bodies, where their `params%new` ran, patch `p11.py`, with `p8.py`
  and `p9.py`): E2, R2, R2FD and MKC identical to the pure move's build. The
  delivered build is identical on the cases the master path does not reach
  (E2S, R2S, R2F, R2F1, MKC, SW1, SW6) and ran E2, R2 and R2FD to the end
  with finite results.
- **Thread counts.** A probe build ran R2 with two parts of six threads and
  a twelve-thread master: the alignment ran at 6, the master was back at 12
  after it, assembly and the commit ran at the master's 12
  (`logs/thread_probe_t11b.log`), finite.
- **Checks after each step:** the Debug build with the fast gate, 15 of 15.
- **Phase-end runs** on the delivered build: SW6 identical; E2S at 24
  threads, MSI with two parts of 12 threads and E2 with two parts of 12
  threads ran to the end with no NaN or infinity (`logs/mt_i11b.log`); the
  full `./compile_debug.sh` with no job running passes, gate 15 of 15 in 16 s
  (`logs/compile_debug_p4_final.log`).
- `lib_stream`, `mini_stream_6vxx` and `mini_stream_1jxy` pass on the
  Debug build (`logs/ctest_stream_p4.log`); the stream pool runs the
  distributed refine2D this phase changed.

Exit: item 11 kept and judged as "The tests" says (the pure move by
identity, the second step on its check build), its evidence above; the
phase-end runs done.

### Phase 5, 2026-10-20: assembly on the master's project, one write per iteration

All three items kept. Files changed: `simple_volassemble_strategy.f90`,
`simple_refine3D_strategy.f90`, `simple_commanders_rec_distr.f90`,
`simple_builder.f90`, `src/main/project/simple_sp_project.f90` and
`simple_sp_project_core.f90`; this plan. Builds `R/scratch/keep/build_i12`
to `build_i14` and `build_c12`; seed patch `p12.py` and the probe
`thread_probe12.py` in `R/scratch/keep/seed_patches`; journal
`R/state/journal_phase_5.md`.

- **Volume assembly on the master's builder** (isolated by ruling). The
  assembly body (`volassemble_exec`) takes the caller's parameters, builder
  and thread count; `run_volassemble` keeps building its own for the
  command-line callers (`reconstruct3D`, the shared-memory refine3D, the
  final reconstruction). The distributed master runs it on its own builder,
  so its particle field receives this iteration's resolution fields in
  memory; the copy of its project that assembly used to hand back is gone,
  and the matching low-pass refresh reads the master's own project.
  Judged on the check build (one `seed_rnd` where assembly's `params%new`
  ran, patch `p12.py`, with the earlier patches) against Phase 4's check
  build: SW1 641 files, CCI 50, PCI 53, CRD 51, SIG 27 identical; MSI
  differs only where the ruling expects: the particle dumps of both steps
  in `cfar` and `lp_est`, and the `LP_ESTIMATED` statistic of both steps,
  which now reports this iteration's nonuniform-filter limit instead of the
  previous one (`R/scratch/compare_cfar.sh` strips the two keys and lists
  every other changed line). The delivered build ran every case it reaches
  to the end with finite results.
- **One write per iteration.** `merge_algndocs` can merge the workers'
  alignment documents into the master's field without writing; assembly and
  the PCG master (the preconditioned conjugate gradient reconstruction) no
  longer write the project on the master's path (a write flag, on by default
  for every other caller); the reprojection model no longer writes, and its
  segment reads happen at the first iteration only. The master writes the
  particle field once before the workers start (with the alignment's sample,
  or as the model left it), again after the merge only when a polish pass
  follows (its workers read the merged field), and once at the end of the
  iteration; the `out` segment keeps its single write before postprocessing
  (which reads it) and does not change after it; the final full write is as
  before. Identity against item 12's build: SW1, MSI 226, CCI 50, PCI 53,
  CRD 51, SIG 27, REC 29, E2 147, R2 32, MKC 19. MSI being identical shows
  that the Nyquist clamp of the resolution in assembly does not bind.
- **The typed entry point and the copied project gone.** `execute_hand_on`
  of the volume-assembly commander and the builder's `spproj_src` option
  (a project copied instead of read) are deleted; every commander's
  `execute(cline)` is now a thin wrapper of a strategy routine. Identity
  against item 13's build on the same cases.
- **Thread counts.** A probe build ran MSI's first step with two parts of six
  threads and a twelve-thread master: assembly ran at the master's 12, the
  alignment at 6, the master was back at 12 after each body
  (`logs/thread_probe_t14.log`), finite.
- **Checks after each item:** the Debug build with the fast gate, 15 of 15.
- **Phase-end runs** on item 14's build: SW6 identical; AD1 passes with its
  check values identical to the base; CR1 passes its 60 checks (its
  distributed stages' values change with the draws, as in Phase 3); FLX
  passes its 95 floors and reaches the weighted branch; E2S at 24 threads,
  MSI with two parts of 12 threads and E2 with two parts of 12 threads ran
  to the end with no NaN or infinity (`logs/mt_i14.log`); the full
  `./compile_debug.sh` on the idle machine passes, gate 15 of 15 in 15 s
  (`logs/compile_debug_p5_final.log`); `lib_stream`, `mini_stream_6vxx` and
  `mini_stream_1jxy` pass (`logs/ctest_stream_p5.log`).

Exit: every item kept and judged as "The tests" says (item 12 on its check
build with only `cfar`, `lp_est` and their statistic differing, items 13
and 14 by identity), its evidence above; the phase-end runs done.

### Phase 6, 2026-10-20: final checks, documents and close

Files changed: `doc/policies/3D/refine3D_policy.md`,
`doc/policies/2D/solve2D_policy.md`,
`doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`,
`doc/policies/3D/reconstruct3D_pcg_policy.md`,
`doc/policies/NU/nonuniform_filtering_policy.md`,
`doc/policies/3D/solve3D_addon_policy.md`,
`doc/implementation_notes/planned/speed_architecture_proposals.md`, the report
`structure2D_and_project_owner_report.md`, and this note, moved from
`planned/` to `completed/`. Journal `R/state/journal_phase_6.md`.

- **The truth suite** (`simple_test_exec test=nightly`, the 19 simulated
  gate entries packed on the 24 cores, the CTest environment) ran on the base
  build and on the final build (item 14's; this phase changes documents only):
  every entry passes on both (base 2,618 s, final 2,720 s), and all 855
  metric rows pass their floors on both
  (`logs/nightly_report_base.tsv`, `logs/nightly_report_final.tsv`, side by
  side with each floor in `logs/nightly_compare.tsv`). The values differ
  only in multi-threaded entries (`cont_refine3D_1jxy`, `flex_pca_blobs`,
  `multistate_assignment`, `multistate_workflows`, `solve3D_addon`, five rows
  of `simulated_workflow_1jxy`), whose draws vary with the threads and, in
  their distributed parts, changed by ruling. No entry fails on either build.
- **Documents.** `refine3D_policy.md`: the strategy routines own probability
  alignment, volume assembly and the sigma2 commit, the commanders are thin
  wrappers, and the distributed master owns its project (its model on its own
  builder, the writes of an iteration). `solve2D_policy.md`: the distributed
  2D master's own project and assembly's class FRCs. The PCG,
  nonuniform-filter, solve3D add-on and multi-state memory notes name the new
  owners; `speed_architecture_proposals.md` records the status.

Exit: the truth suite on both builds with every entry passing, the documents
written, the report beside this note, the note in `completed/`.
