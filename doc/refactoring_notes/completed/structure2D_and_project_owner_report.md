# The 2D structure aligned with 3D, and one project owner on the master: report

Report of the unattended run `structure2D_owner` (2026-10-20, on the Oracle
Linux test machine). The plan, with every number, log path and decision in its
"Progress" section, is `structure2D_and_project_owner.md` beside this report.
The diffs of the phases are in `~/agent_runs/structure2D_owner/review/phase_N.diff`
(N = 0 to 6); `R` below is `~/agent_runs/structure2D_owner`.

## What the run learned first

- The last three speed commits are byte-identical to their parents on this
  machine's x86 build too, as they were on the Mac.
- Every step the restructuring touched could be kept byte-identical at one
  thread, except where a ruling accepted a change: the moved steps no longer
  reseed the random generator (their check builds with the reseeds put back
  are byte-identical), and `cfar` and `lp_est` reach the master one iteration
  earlier.
- Two defects of the base were found: a fresh 2D start with fractional update
  and `fillin=yes` never ended (fixed), and shared-memory and distributed 2D
  restoration use different class FRCs (Fourier ring correlations) for the
  even/odd low-resolution averaging (kept as it is, for the maintainer).

## What changed, phase by phase

**Phase 0, the base and the x86 check** (`phase_0.diff`). Release and Debug
builds; the identity cases made on the base build and run twice (all
deterministic at one thread); the exclusions declared. The x86 check:
`aea5b1846` against `508d99066`, `af13bc828` against `aea5b1846`,
`f078cf683` (`cache=no`) against `0c9542e50` (`cache=yes`): all identical,
the probability-table part files compared through lossless decoders across
their old and new layouts. One code fix: `refine2D` started fresh
(`startit=1`) with `update_frac` below 1 and `fillin=yes` reused its first
particle subset for ever while `fillin` waited for every particle to be
assigned, so it never ended; a fresh start with `fillin` now samples by
update count (`simple_matcher_smpl_and_lplims.f90`).

**Phase 1, the 2D structure aligned with 3D** (`phase_1.diff`). The 2D matcher
(`refine2D_exec`) searches, accumulates the class sums and writes the
orientations, and a distributed worker its contribution; the strategy restores
the class averages in shared memory (`restore_cavgs_inmem`) and one routine
decides convergence in both modes (`check_refine2D_convergence`). `refine2D`
is a public `simple_exec` program (required keys `projfile` and `mskdiam`);
the workers stay on `simple_private_exec`. `solve2D_policy.md` describes it.

**Phase 2, the 3D bodies out of the commanders** (`phase_2.diff`). Volume
assembly moved to `simple_volassemble_strategy.f90`, `prob_align` and
`prob_align_neigh` to `simple_prob_align_strategy.f90`, the sigma2 commit to
`simple_strategy_sigma2.f90`, unchanged; assembly returns the names of the
volumes it wrote instead of writing `params%vols`.

**Phase 3, the 3D master owns its project** (`phase_3.diff`). The distributed
refine3D master makes its reprojection model (the references the workers
search) on its own builder instead of a second kept builder, and no longer
reads the whole project after it; `prob_align` and the sigma2 commit run on
its `(params, build)` with explicit thread counts; the `execute_from_project`
entry points are gone.

**Phase 4, the 2D master owns its project** (`phase_4.diff`). `prob_align2D`
and class-average assembly became strategy routines (the latter in the new
`simple_cavgassemble_strategy.f90`), then the distributed refine2D master ran
them and the sigma2 commit on its own `(params, build)`; assembly writes only
the segments it changes.

**Phase 5, assembly on the master's project, one write per iteration**
(`phase_5.diff`). Volume assembly runs on the distributed refine3D master's
builder; the master writes its particle field once before the workers and once
per iteration (the merge of the workers' documents in memory, no project writes
from assembly or the PCG master on that path, the rereads gone);
`execute_hand_on` and the builder's `spproj_src` option are deleted. Every
commander involved is now a thin wrapper of a strategy routine.

**Phase 6, final checks and documents** (`phase_6.diff`). The truth suite on
the base and the final build; `refine3D_policy.md`, `solve2D_policy.md`, the
PCG, nonuniform-filter, solve3D add-on and multi-state memory notes; the status
of `speed_architecture_proposals.md`; this report; the plan moved to
`completed/`.

## Decisions taken under delegation

Each is recorded in the plan's progress section of its phase.

- MSI, the three-state identity case, starts from E2's 1,500 simulated
  particles with their simulated poses: on x86 the fixture's 82 particles left
  a state with one particle and the case's second step refused to start
  (Phase 0; the plan's table of cases says so).
- R2F and R2FD continue from E2's last stage (`startit=2` with its carried
  class sums), because a fresh start never blends carry-over; R2F1 keeps the
  fresh-start form as a check of the sampler fix (Phase 0).
- The graphical interface's metadata (`metadata.json`) of SW1's extraction
  step is written by the master and its worker at the same moment; the
  comparison skips that file where two processes wrote previews (Phase 0).
- The public `refine2D` versus `refine2D_distr` comparison maps the execution
  directory in its relative form too (`1_refine2D/` against `run/`, in the star
  files) and checks the untouched input project separately; no file excluded
  (Phase 1, `R/scratch/exclusions.tsv`).
- Volume assembly's command-line entry moved with its body, so that the
  strategy modules stop importing the commander; the commander is the thin
  wrapper (Phase 2).
- The probability-alignment bodies take the in-process table commander as an
  optional argument of the abstract commander type, because importing it from
  the module that calls the bodies would be circular (Phase 2).
- The master's model reads the segments the kept builder read, before the
  model, and names its own project (`update_projinfo`); from Phase 5 only at
  the first iteration (Phases 3 and 5).
- Shared-memory refine2D and refine3D stay on the commanders' `execute(cline)`
  (Phases 3 and 4); the command-line class-average assembly keeps its full
  project write, which keeps the file's `projinfo` current (Phase 4).
- Class-average assembly starts its class FRCs afresh, as a new toolbox did,
  to stay identical (Phase 4).
- The `out` segment keeps its single write before postprocessing, which reads
  it from the file, instead of a second write at the end of the iteration
  (Phase 5).
- The matching low-pass refresh takes a flag to read the master's own project
  instead of a copied project argument (Phase 5).

## Deviations from the plan and incidents

- **The fixture was regenerated** (Phase 3). The first build of the master's
  own builder left out `update_projinfo`; a run started from a project in
  another directory names that project, and the master's writes went to the
  fixture's project. The identity check caught it; the fixture was
  regenerated on the base build at its original path and verified (every
  file but that project equal to the damaged copy, the project equal to SW1's
  reference run in every segment). The plan says the fixture is never rerun.
- **FLX is not deterministic** (Phase 2): it runs its per-state refinements at
  two threads whatever `nthr` is, so it is judged by its floors, not its
  values. All 95 pass on every build.
- **REC's reference was made again** on Phase 4's build in Phase 5, because the
  older REC trees had lost their large files in clean-ups.
- **Thread counts** were shown with scratch-only probe builds (log lines of the
  team sizes), because the code logs none.

## The evidence against the base

- Identity: every case each item reaches, at one thread with a fixed seed,
  byte for byte against the build before the item (project files through
  their dumps, table part files through decoders); the seed items on check
  builds whose reseeds are put back by scratch-only patches
  (`R/scratch/keep/seed_patches/p8.py`, `p9.py`, `p11.py`, `p12.py`).
- After every item: the Debug build with the fast gate (15 of 15). At every
  phase end: SW6, the full Debug build on an idle machine, `lib_stream` and
  the two `mini_stream` entries (the stream's 2D pool runs the distributed
  refine2D), multi-threaded runs to the end; AD1 and CR1 at the ends of Phases
  1, 3 and 5 (AD1's check values identical to the base throughout; CR1's
  distributed stages change with the draws from Phase 3 on, all checks pass).
- The truth suite (`simple_test_exec test=nightly`, 19 entries) on the base and
  the final build: every entry passes on both, and all 855 metric rows pass
  their floors on both (`R/scratch/logs/nightly_compare.tsv`); the values
  differ only in multi-threaded entries, whose draws vary with the threads and,
  in their distributed parts, changed by ruling.

## Open for the maintainer

- **The FRC difference between the modes**: shared-memory refine2D reads the
  previous iteration's class FRCs from file for the even/odd low-resolution
  averaging; distributed assembly starts from zero FRCs at every assembly.
- **The metadata race**: the extraction master and its worker both write the
  graphical interface's `metadata.json` (`src/utils/gui`).
- **The graphical interface** (`nice/`, not edited): class selection is
  accepted only from `solve2D` (`nice/nice_lite/views/job_builder_views.py`),
  so `refine2D` results cannot be selected there yet. `doc/user_guide/user_guide.md`
  does not list `refine2D`.
- **Stream code**: nothing needed; the stream pool runs the distributed
  refine2D through its commander and all stream tests pass.
- **`refine2D_nano`**: with `prg=refine2D` registered, an in-process call
  honours `mkdir` and project discovery; `refine2D_nano` sets both, but a user
  `mkdir=no` without a project file would now look for one.
- **The `simple-refine3d` skill** says volume assembly is intentionally
  executed in the commander. Proposed text: "Volume assembly is a strategy
  routine (`volassemble_exec` in `simple_volassemble_strategy.f90`) on a
  caller's `(params, build)`; the distributed refine3D master runs it on its
  own builder and owns the project, writing the particle field once before the
  workers and once per iteration; `commander_volassemble` is a thin wrapper for
  command-line callers. `prob_align`, `prob_align_neigh` and the sigma2 commit
  follow the same pattern." The `simple-solve2d` skill could add that the 2D
  strategy restores the class averages and decides convergence, and that
  `refine2D` is public.

## Commit message

    The 2D structure aligned with 3D; one project owner on the master

    The 2D matcher searches and accumulates; the 2D strategy restores the class
    averages and decides convergence; refine2D is a public program. The bodies
    of prob_align, prob_align_neigh, prob_align2D, volume assembly,
    class-average assembly and the sigma2 commit are strategy routines on the
    distributed master's (params, build); the master owns its project and writes
    it once per iteration. Fixes a fresh 2D start with fillin that never ended.
