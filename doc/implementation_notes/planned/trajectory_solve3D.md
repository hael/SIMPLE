# `trajectory_solve3D`: the nanoparticle pipeline in one program, with a second gate on the class averages

An unattended run on the Dell (`trajectory_solve3D`), 2026-10-29. It builds the one-command workflow from the
graphene-subtracted frames of a tracked particle to its time-integrated map, adds a second rejection of class
averages by their agreement with a start model, sets that gate on the maintainer's labelled trajectories, and
runs every platinum trajectory of `/data2/NanoX_trajs` through it beside the products of the previous run
(`nano_pipeline_pt`, 2026-10-29, which ran the same stages by hand and found the rejection leaving a tail of
poor class averages in most trajectories).

## 1. What the maintainer fixed (2026-10-29)

- One program does it all: `trajectory_solve3D`, a new SINGLE program (`single_exec`), named for what it does,
  a de novo map from a trajectory; `solve3D_nano`, the thin wrapper around `solve3D`, stays as it is.
- The second gate aligns the class averages to the start model simulated from the diameter estimate (the FCC
  particle `analysis2D_nano` already simulates), not to a map made from the frames: a class average that
  matches no projection of a crystalline particle of the measured size at lattice resolution is rejected, and
  nothing in that judgement comes from the data being judged.
- The gate is a relative cut on the correlation, set once on the three labelled trajectories for recall
  first (as the lattice gate was: for these time series losing good frames costs more than carrying some junk),
  then fixed as a constant, with the same kind of override key the lattice gate has.
- The stages otherwise as the previous run ran them: flex denoising with its defaults, `refine2D_nano`, the
  lattice gate (`quality_model=nano`), `center2D_nano`, `autorefine3D_nano` without the atomicity prior and
  with `icm=no`.
- Bugs are fixed as they are found, in the phase that finds them. The maps are not judged by the run beyond
  the numbers it tables; the maintainer validates them.

Terms. A trajectory is one tracked particle: a stack of frames, graphene-subtracted (`subtr_frames.mrc`), 0.358 Å.
`refine2D_nano` groups the frames into class averages that are time chunks of about twenty frames. The lattice
gate accepts a class average whose strongest Bragg peak stands at least `10^0.65` times above its ring. The
start model is the FCC platinum particle `simulate_nanoparticle` builds for the estimated diameter. Aligning a
class average to a model means finding the model projection it matches best (`refine3D_nano` on the `cls3D`
field, one iteration) and its correlation to it in the band 3 to 1 Å (the program's `hp` and `lp`). The time-
integrated map is what `autorefine3D_nano` refines from the kept frames. The labelled trajectories are
71-17/traj04, 71-6/traj01 and 71-6/traj04, whose class averages the maintainer sorted by eye in October 2026
(`cls2D_field.txt`, `ptcl2D_field.txt` beside the frames; the previous run on the rejection used them).

## 2. The program

`single_exec prg=trajectory_solve3D projfile=<p> stk=<subtr_frames.mrc> smpd=0.358 nthr=24`, on a project made
by `new_project`. Optional: `mskdiam` and `moldiam` (else estimated), `lp` for the 3D (default 1.5), `thres`
for the lattice gate, the second gate's override key, `pca_mode`, `halfwin`, `neigs`, `eo_pool` for the
denoiser. The program runs, in its execution directory, each stage in its own numbered subdirectory with its own
project copy, as the SINGLE programs do, so that every intermediate is where the previous run left it:

1. **Denoise and import.** `trajectory_denoise` on `stk` (flex defaults), `trajectory_import` of the denoised
   stack into the project.
2. **Diameter, first estimate.** Ten time-chunk averages (`make_cavgs ncls=10 tseries=yes nparts=1`),
   `estimate_diam` (`lp=3`, `mskdiam=0`): `mskdiam = 1.5 x max_diam`, `moldiam = min_diam`. Where that mask
   exceeds the box, the parameter layer's box default applies, as now (recorded in the log).
3. **2D.** `refine2D_nano` with that `mskdiam` and its defaults.
4. **Gate 1, the lattice gate.** `model_cavgs_rejection quality_model=nano mskdiam=<mskdiam>` without `prune`:
   the class states and the quality table; the accepted and rejected montages it writes.
5. **Start model.** `simulate_nanoparticle element=Pt moldiam=<moldiam> box=<box>`: `startvol.mrc`.
6. **Gate 2, the model gate.** The class averages that passed gate 1 are aligned to `startvol.mrc`
   (`refine3D_nano oritype=cls3D maxits=1 lp=1.0 trs=10 mskdiam=<mskdiam>`, as `autorefine3D_nano`'s own final
   pass does); the correlation of each to its best reprojection is read from `cls3D`. The cut: a class is
   rejected when its correlation lies below `median - K_MODEL x 1.4826 x MAD` of the gate-1 survivors' correlations
   (median and median absolute deviation over those classes; `K_MODEL` the constant phase 0 sets). States go to
   `cls2D`, `cls3D`, `ptcl2D` and `ptcl3D` through `map_cavgs_selection`; the project is pruned
   (`prune_particles`). Products: the per-class table (class, population, lattice peak, gate-1 state,
   correlation, best orientation, gate-2 state), and the triplets class average / reprojection for the classes
   gate 2 rejected and for a sample of those it kept, for the maintainer's eye.
7. **Centering.** `center2D_nano ncls=10` on the pruned project.
8. **Diameter, second estimate**, on the centred chunk averages; it sets `mskdiam` and `moldiam` for 3D (the
   start volume is simulated again if `moldiam` changed by more than 5 percent).
9. **3D.** `autorefine3D_nano vol1=startvol.mrc pgrp=c1 mskdiam=<mskdiam> lp=1.5 icm=no`, its defaults
   otherwise, no `atm_prior`. Its products as now (`final_results/`, the iterations' text products, the class
   averages aligned to the final map with their triplets, the per-frame correlations).

The workflow writes `trajectory_solve3D.txt` in its directory: one line per stage with its wall time and the
numbers of section 4, and a last line with the map's path. A stage that fails stops the program with the
stage's name in the error; the stages already done are kept so that a rerun resumes.

Ownership: the commander beside `autorefine3D_nano` in `single_commanders_nano3D.f90`, driving the existing
commanders (`single_commanders_trajectory`, `single_commanders_nano2D`, `simple_commanders_cavgs`,
`simple_commanders_sim`, `simple_commanders_imgproc`) as `analysis2D_nano` and `autorefine3D_nano` do; the gate-2
rule (the relative cut and its constant) in `src/main/cavg_quality/` beside the lattice gate, with its tester
case, and the policy note `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` describing the second
gate; the UI entry in `single_ui_trajectory.f90` with every key named above, the router case in
`single_exec_trajectory.f90`, the parameter validation in `simple_parameters_phases.f90`. Follow
`doc/policies/layering_and_project_ownership_policy.md`: the workflow reads and writes the project through its
API, one owner per stage.

## 3. Setting the gate on the labelled trajectories (phase 0)

The labels are classes of an earlier 2D run (scattered classes of an earlier denoiser), not of this
workflow's time chunks, so the gate is set on the labelled objects themselves and checked in situ afterwards.

1. For each labelled trajectory: flex-denoise `subtr_frames.mrc` (the stack of the previous run is
   bit-identical and can be regenerated in 80 s); rebuild the maintainer's classes on that stack from the
   stored assignment and alignment (`import_particles stk=<flex stack> smpd=0.358 ctf=no
   oritab=ptcl2D_field.txt`, `make_cavgs ncls=<classes> refs=cavgs_labelled.mrc`), and take their labels from
   `cls2D_field.txt` (`state`). Check the populations against the delivered ones class by class, as the previous
   rejection run did.
2. Gate 1 on those class averages (`model_cavgs_rejection quality_model=nano`, no prune), the start model from
   `estimate_diam` on them, gate 2's alignment and correlations as in section 2.
3. The tables, per labelled set and pooled: for every class its label, population, lattice peak, gate-1 state,
   correlation. Then the three gates against the labels: the lattice gate alone; the correlation cut alone at
   `K` = 1.5, 2, 2.5, 3, 4; the lattice gate followed by the cut. For each: recall of the maintainer's selected
   classes, the fraction of his rejected classes removed, balanced accuracy, and the area under the curve of the
   correlation alone and of the lattice peak alone (the previous run measured the latter at 0.97 to 0.99).
4. The choice: the smallest `K` in the list (the strictest cut) at which, with the lattice gate in front, the
   recall of the maintainer's selected classes is at least 0.98 on every labelled set; if no `K` in the list
   reaches it, `K` = 4 and the fact reported. `K_MODEL` is that constant, recorded with the table in the
   Progress entry and in the policy note. It is set once here; it is not revisited on the downstream numbers.
5. In situ, on the time-chunk classes of the workflow (phase 2): each class gets the majority label of its
   frames (from `ptcl2D_field.txt`'s per-frame state through `indstk`) and a purity; the gates' recall and
   removal on the classes with purity 0.8 or more are reported beside the labelled-set numbers. A report, not a
   second setting of `K`.

## 4. What is kept, and the comparison

Per trajectory, under `scratch/keep/pipeline/<dataset>/<trajectory>/`: the workflow's directory less the bulky
intermediates (the denoised stack of boxes above 120 dropped after the run, as before; the per-iteration
volumes dropped; every project, log, table, montage and triplet stack kept; `final_results/` kept whole).

`scratch/keep/pipeline/table.tsv`, one row per trajectory, the columns of the previous run's table (frames,
denoise numbers, 2D, gate 1, diameters, centering, 3D wall time, resolution, atoms, `detect_atoms` fit) plus:
classes gate 2 rejected and their frames; the correlation median, spread and cut; frames kept after both gates
and the fraction; and, beside each 3D number, the previous run's value for the same trajectory from
`~/agent_runs/nano_pipeline_pt/scratch/keep/pipeline/table.tsv` (read only). From the 3D's per-frame
correlations (`ptcls_vs_reprojs_corrs.csv`): the fraction of frames below `median - 3 x 1.4826 x MAD` and
whether those frames cluster in time (the longest run of consecutive such frames), reported as the measure of
the remaining tail, not acted on.

## 5. Rules of the run

- Fix as you go, as the previous run: a stage that fails or misbehaves is fixed in its owning module (section
  7; a file not listed is added with the reason before it is edited), with a check in the area's tester where
  one exists and the fast gate rerun; the trajectories already past that stage are rerun from it.
- The constant `K_MODEL` is set in phase 0 on the labelled sets and not touched afterwards; no threshold or
  default of the pipeline moves to let a trajectory through.
- One trajectory at a time at 24 threads, smallest first (the previous run's inventory order). The program is
  resumable per stage.
- The tests: `trajectory_solve3D` is exercised in `single_workflow` (no new CTest entry): the `fcc` suite's
  simulated trajectory, with three contiguous windows of twenty frames replaced by pure noise before the
  denoiser, goes through `trajectory_solve3D` instead of the hand chain, and the test asserts that at least
  90 percent of the noise frames are gone after the gates, that at most 5 percent of the true frames are, and
  its present floors on the final map; the `wurtzite` suite keeps the hand chain, so both forms stay covered.
  The gate-2 rule has a closed-form case in the cavg-quality tester (a given correlation list, the expected
  cut). The density of the fixture is not changed. `check_test_registry.py` and `check_descr.py` whenever a
  table, router or header changes; `SIMPLE_CTEST_BUDGET` unchanged.
- Repository rules as always: no backwards compatibility; `THROW_HARD` arguments on one line; no goto; no
  short-circuit on allocatable elements; dead code met is deleted; new OpenMP code at 11 and 24 threads; mode
  changes checked; scripts of the run under `scratch/keep/`, not in the repository; no commits.

## 6. Phases

**Phase 0, the gate on the labelled sets.** Build, fast gate. The labelled trajectories' files present
(`cls2D_field.txt`, `ptcl2D_field.txt` beside `subtr_frames.mrc`; if absent, stop and say so). The study of
section 3, steps 1 to 4, with pure-Python tables; `K_MODEL` chosen and recorded. Ends DONE with the tables and
the constant.

**Phase 1, the program.** `trajectory_solve3D` as section 2 describes, its UI, router, parameters, the gate-2
rule and its tester case, the policy note, the `single_workflow` change of section 5; the fast gate, both
`single_workflow` entries, `unit_single`. Then one real trajectory (71-17/traj01) end to end through the
program, its products compared stage by stage with the previous run's (the same denoised stack bit for bit,
the same 2D within the run-to-run variation it had, the gate-1 states identical). Ends DONE when the test
passes and the trajectory has its map.

**Phase 2, the ten trajectories.** All of `/data2/NanoX_trajs` through the program, the table of section 4
after each, the in-situ check of section 3 step 5 on the three labelled ones. Fix as you go. Ends DONE with
every trajectory mapped or stopped with its diagnosis.

**Phase 3, the report.** `doc/implementation_notes/planned/trajectory_solve3D_report.md`, reading on its own:
what the program does, how the gate was set (the labelled-set tables), what gate 2 removed on the ten
trajectories and what that did to the 3D numbers beside the previous run's, the remaining tail per frame, every
fix, where the products are, what is left open.

## 7. File table

A phase edits only its files; a file the work turns out to need is added here, with its phase and a one-line
reason, before it is edited. Testers (`*_tester.f90`) of a listed Fortran file are covered by its row.

| File | What changes | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/trajectory_solve3D.md` | this note: Progress, file-table rows | each |
| `doc/implementation_notes/planned/trajectory_solve3D_report.md` | the report | 3 |
| `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` | the second gate and its constant | 1 2 |
| `src/main/commanders/single/single_commanders_nano3D.f90` | the `trajectory_solve3D` commander; fixes | 1 2 |
| `src/main/ui/single/single_ui_trajectory.f90` | the UI entry | 1 2 |
| `src/main/exec/single_exec_trajectory.f90` | the router case | 1 2 |
| `src/main/params/simple_parameters.f90` | the gate-2 override key | 1 2 |
| `src/main/params/simple_parameters_parse.f90` | the same | 1 2 |
| `src/main/params/simple_parameters_phases.f90` | its validation and the program's | 1 2 |
| `src/main/cavg_quality/` | the gate-2 rule and its constant, with the tester case | 1 2 |
| `src/main/commanders/test/simple_commanders_test_single.f90` | `single_workflow`: the `fcc` suite through the program with noise windows | 1 2 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the `single_workflow` description | 1 |
| `production/CMakeLists.txt` | the `single_workflow_fcc` comment | 1 |
| `doc/policies/test_environment_policy.md` | the `single_workflow` claim | 1 |
| `src/main/commanders/single/single_commanders_trajectory.f90` | trajectory_denoise, trajectory_import: fixes | each |
| `src/main/flex/simple_flex_traj_denoise.f90` | the denoiser: fixes | each |
| `src/main/commanders/single/single_commanders_nano2D.f90` | center2D_nano, refine2D_nano: fixes | each |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | model_cavgs_rejection: fixes | each |
| `src/main/strategies/search/simple_strategy2D_tseries.f90` | the time-series 2D search: fixes | each |
| `src/main/strategies/search/simple_strategy2D_alloc.f90` | the fresh-start rule: fixes | each |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | refine2D as driven here: fixes | each |
| `src/main/strategies/parallelization/simple_make_cavgs_strategy.f90` | the chunk averages: fixes | each |
| `src/main/nano/simple_nanoparticle.f90` | detect_atoms inside the loop: fixes | each |
| `src/main/nano/simple_nano_density.f90` | the density model: fixes | each |
| `src/main/commanders/simple/simple_commanders_atoms.f90` | detect_atoms commander: fixes | each |
| `src/main/commanders/simple/simple_commanders_sim.f90` | simulate_nanoparticle: fixes | each |
| `src/main/commanders/simple/simple_commanders_imgproc.f90` | estimate_diam: fixes | each |
| `src/main/project/simple_sp_project_ptcl.f90` | map_cavgs_selection, prune: fixes | each |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | model_cavgs_rejection records each nano class's raw lattice peak in `cls2D`, for the program's class table | 1 |
| `src/main/params/simple_parameters_registry.f90` | ruling 4: a refused file format stops with `THROW_HARD` (it exited with status 0) | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | ruling 4: `refine3D oritype=cls3D` stops with a message when `cls3D` does not hold one entry per class average | 1 |

## 8. Progress

(one entry per phase, written by the run)

### Phase 0, the gate on the labelled sets (2026-10-29, DONE)

**In short.** The second gate's constant is set: `K_MODEL = 3`. With the lattice gate in front, 3 is the
strictest cut in the plan's list that still keeps at least 98 percent of the maintainer's selected classes on
each of the three labelled trajectories. On these sets the gate removes little at that constant: nothing on
71-17/traj04 or 71-6/traj04, and 9 classes on 71-6/traj01, all of them classes the maintainer rejected that sit
just above the lattice threshold. The correlation to the start model separates good from bad classes as well
as the lattice peak does (area under the curve, AUC, 0.96 to 0.99 against 0.975 to 0.99). But a cut placed
relative to the median and spread of the survivors stays loose, because the survivors still hold a second,
junk population that widens the spread. One bug was found and fixed: `make_cavgs` did not write the
class-average orientation field (`cls3D`) it creates, so the alignment against the start model crashed on its
project. Ruling 2 (the second gate, set for recall first on the labelled sets) is applied here for the first
time, as written.

**Build.** Release build with tests (`./compile_clean.sh`), the same toolchain as the maintainer's
`~/src/SIMPLE/build` and the previous run: gfortran 15.2.1 from gcc-toolset-15, system FFTW, CMake 3.26.5.
Every job sources `scratch/keep/env.sh`. The fast gate (the 15 quick unit tests run at build time) passed 15 of
15 after the clean build and again after the fix (`scratch/keep/build_clean.log`,
`scratch/keep/fast_gate_after_fix.log`). Paths in this entry are under
`/home/elmlundho/agent_runs/trajectory_solve3D/`.

**The labelled files and the rebuild (section 3, step 1).** `cls2D_field.txt` and `ptcl2D_field.txt` are present
beside `subtr_frames.mrc` for all three trajectories. The flex denoiser with its defaults took 80, 74 and 60 s.
Its parameter table equals the previous run's on all three, and the denoised stacks of 71-6/traj01 and
71-6/traj04 are bit-identical to the previous run's (that of 71-17/traj04 had been deleted). The classes were
rebuilt on these stacks from the stored assignment and alignment (`scratch/keep/gate/run_set.sh`):

1. `simple_exec prg=new_project projname=lab`
2. `simple_exec prg=import_particles projfile=lab.simple stk=<flex stack> smpd=0.358 ctf=no kv=300 cs=2.7 fraca=0.1 oritab=ptcl2D_field.txt`.
   The program requires `kv`, `cs` and `fraca` even with `ctf=no`.
3. `simple_exec prg=make_cavgs projfile=1_import_particles/lab.simple ncls=<649|684> refs=cavgs_labelled.mrc nparts=1 nthr=24`.

The populations equal the delivered ones class by class on all three sets (649, 684 and 684 classes, no
mismatch). Each frame's class, in-plane angle, shift, stack index and even/odd half equal the delivered table's.
The frame count equals the stack's (12,971, 13,662 and 13,662). The import sets every frame's state to 1, so
the maintainer's labels never enter the project; they are read from `cls2D_field.txt` by the table script. The
empty classes are 25 and 122 (71-6/traj01) and 9 and 45 (71-6/traj04). They are left out of every number below:
the lattice gate rejects them outright, and their correlation to a projection is undefined.

**The gates (section 3, step 2), as run on each set.**

- Diameter: `single_exec prg=estimate_diam stk=cavgs_labelled.mrc smpd=0.358 mskdiam=0 lp=3`, with
  `mskdiam = 1.5 x max` and `moldiam = min`.
- Gate 1: `simple_exec prg=model_cavgs_rejection quality_model=nano mskdiam=<mskdiam>`, without prune.
- Start model: `single_exec prg=simulate_nanoparticle element=Pt moldiam=<min> box=<box>`.
- Gate 2's alignment: `single_exec prg=refine3D_nano oritype=cls3D maxits=1 lp=1.0 trs=10 icm=no volrec=no vol1=startvol.mrc pgrp=c1 mskdiam=<mskdiam>`,
  on a copy of the rebuilt project. Every non-empty class is aligned, so that the cut can also be tabled
  without the lattice gate in front; the empty classes are first set to state 0
  (`simple_exec prg=selection oritype=cls3D deselfile=`). One pass takes 6 to 10 s.

What the alignment does was read from the code before the numbers were trusted. With `refine3D_nano`'s
`refine=neigh`, the first and only iteration is a global search: 500 coarse directions, then refinement of the
best on the 2,500-direction grid, then in-plane angle and shift within `trs`. A class's `corr` in `cls3D` is its
best correlation (`objfun=cc`) to a reprojection of the start model in the band 3 to 1 Å. The model is
soft-masked at `mskdiam` and low-passed to 1 Å. `volrec=no` skips a reconstruction from the class averages,
which does not change the correlations. The reprojections for the montages come from
`simple_exec prg=reproject ... oritab=<cls3D orientations> mskdiam=<mskdiam> pgrp=c1 nspace=2500`. The user
interface (UI) requires `nspace` even when the orientations are given.

**A decision on the start model.** The diameter estimate over the 649 to 684 scattered labelled classes has
wider extremes than the program's own estimate over ten time-chunk averages of the same stack. Labelled
classes, min/max: 20.9/38.6, 16.6/29.4 and 16.7/32.6 Å. Ten chunks: 24.0/27.4, 19.3/23.9 and 17.7/28.9 Å,
identical to the previous run's. The plan says the start model comes from the estimate on the labelled
classes, and that is the primary model. The chunk-based model, the one the program will build, was run beside
it as a check, with the rule that if the two chose different constants the larger (more permissive) one would
be taken. Both chose 3. Where `1.5 x max` exceeds the box (four of the six runs), the parameter layer falls back
to its box default and warns, as section 2 provides.

**Fix: `make_cavgs` writes the class-average fields it makes congruent** (`simple_make_cavgs_strategy.f90`).
Exposed by the first alignment on 71-17/traj04: `refine3D_nano oritype=cls3D` stopped with a segmentation fault
(in `get_state` under `sample4update_all`, from the backtrace through `addr2line`). Registering the class
averages (`add_cavgs2os_out`) makes `cls2D` (its stack index) and `cls3D` congruent with them in memory. The
shared-memory bookkeeping of `make_cavgs`, however, wrote `cls2D` before that step and never wrote `cls3D`. A
project made by `make_cavgs` alone therefore registered class averages in `out` with no `cls3D` on disk, and
`refine3D` read past the empty field. The fix writes `cls2D`, `cls3D` and `out` after the registration. The
distributed path assembles through `cavgassemble_exec` and is unchanged. `refine2D_nano` projects, which the
program's gate 2 will read, already carry `cls3D`, so the program itself was not affected. No tester exists for
the `make_cavgs` strategy. The checks are the rerun (the alignment ran on all three sets), the fast gate,
`unit_project`, and both SINGLE workflow tests, whose `center2D_nano` drives `make_cavgs`: `single_workflow_fcc` passed in 660 s and `single_workflow_wurtzite` in 438 s (`scratch/keep/ctest_single_workflow.log`). Left
open, outside this phase's files: `refine3D` on `oritype=cls3D` has no check that `cls3D` holds as many entries
as there are class averages, so a malformed project still crashes rather than stopping with a message.

**The tables (section 3, step 3)** (`scratch/keep/gate/`, by `tables.py`).

- Per set `<set>/classes.tsv`: class, label, population, lattice peak, gate-1 state, correlation, orientation and
  shift.
- Per set `<set>/gates.tsv`: for each gate and K, the recall of the selected classes, the fraction of the
  rejected classes removed, the balanced accuracy, the classes kept, the frame-weighted recall, and the cut's
  median, robust spread and value.
- Pooled: `gates_pooled.tsv`, and the AUCs in `auc.tsv`.
- The chunk-model check: `<set>/model_chk/` and `sensitivity_chk/`.

Recall is the fraction of the maintainer's selected classes a gate keeps. Removal is the fraction of the classes
he rejected that it removes. Balanced accuracy is their mean.

| Set | classes (selected / rejected) | lattice gate: recall / removal | + cut at K = 3: recall / removal | gate 2 removes | AUC corr / lattice peak |
| --- | --- | --- | --- | --- | --- |
| 71-17/traj04 | 649 (185 / 464) | 0.984 / 0.735 | 0.984 / 0.735 | 0 | 0.990 / 0.990 |
| 71-6/traj01 | 682 (425 / 257) | 0.988 / 0.444 | 0.988 / 0.479 | 9 (all rejected by him) | 0.959 / 0.975 |
| 71-6/traj04 | 682 (264 / 418) | 0.985 / 0.670 | 0.985 / 0.670 | 0 | 0.982 / 0.983 |

Pooled over the three sets, the lattice gate followed by the cut, by K (recall / removal / balanced accuracy,
then the lowest recall over the sets):

| K | 1.5 | 2 | 2.5 | 3 | 4 | lattice gate alone |
| --- | --- | --- | --- | --- | --- | --- |
| recall / removal / BA | 0.929 / 0.786 / 0.857 | 0.950 / 0.725 / 0.837 | 0.982 / 0.684 / 0.833 | 0.986 / 0.653 / 0.820 | 0.986 / 0.645 / 0.816 | 0.986 / 0.645 / 0.816 |
| lowest set recall | 0.871 | 0.913 | 0.979 | 0.984 | 0.984 | 0.984 |

The correlation cut alone, without the lattice gate in front, keeps every selected class and removes at most
1.5 percent of the rejected ones at any K. The rejected classes are the majority on two sets, so the median sits
in the junk population and the cut falls below nearly everything. The pooled AUC of the correlation is 0.982
(lattice peak 0.984). Among the gate-1 survivors it is 0.998, 0.944 and 0.969 (pooled 0.970). The 9 classes
gate 2 removes on 71-6/traj01 have lattice peaks of 0.65 to 0.82, just above the lattice threshold of 0.65, and
correlations of 0.19 to 0.23 against the survivors' median of 0.71.

**The choice (section 3, step 4).** `K_MODEL = 3`. The rule takes the smallest K in the list at which, with the
lattice gate in front, the recall of the selected classes is at least 0.98 on every labelled set. K = 2.5 misses
on 71-6/traj01 (0.979; 0.965 with the chunk model). K = 3 reaches 0.984, 0.988 and 0.985, the lattice gate's
own recall: the cut adds no loss of selected classes on any set. The chunk-model check gives the same choice:
there K = 3 removes 25 classes on 71-6/traj01, 2 of which the maintainer selected, for a recall of 0.984. At
this constant the second gate is a safety net for the junk just above the lattice threshold, not a large second
filter. That follows from the rule's form, a median-relative cut over a population that is still bimodal after
the lattice gate. The form and the recall-first setting are the maintainer's, fixed by ruling 2, and are not
revisited here. Phase 1 writes the constant into
`doc/policies/sieving_and_rejection/model_cavgs_rejection.md` with the rule.

**Montages for the maintainer** (`scratch/keep/gate/<set>/montages/`). At K = 3 they show, as normalised
class-average/reprojection pairs, the classes the lattice gate rejects (`lattice_rejected_pairs.mrc`), the classes
the cut removes behind it (`lattice+cut_K3_rejected_pairs.mrc`, 71-6/traj01 only, since it is empty elsewhere),
and the classes the cut alone would remove (`cut_alone_K3_rejected_pairs.mrc`, 4 classes on 71-17/traj04). Each
stack comes with a `.tsv` naming every pair's class, label, population, lattice peak, correlation and shift. As
in `autorefine3D_nano`'s triplets, the reprojection is not shifted by the class's shift.

**Kept and deleted.** Kept under `scratch/keep/gate/`:

- the scripts (`run_set.sh`, `run_all.sh`, `denoise_all.sh`, `tables.py`, `montages.py`, `mrcpairs.py`);
- every log, project file and table;
- the rebuilt class averages, the start model and its reprojections (the primary model only);
- the montages.

The denoised stacks (regenerated bit for bit in about 80 s) and the other bulky intermediates were deleted after
their numbers were taken (`scratch/keep/gate/deleted_bulky.txt`).

**Files changed in this phase.** `src/main/strategies/parallelization/simple_make_cavgs_strategy.f90` (the fix) and
this note (this entry). `git diff --summary` shows no mode changes. No suite table, tester registration, router or
file header changed; `check_test_registry.py` and `check_descr.py` pass all the same.

**Exit items.** Build and fast gate: done, 15/15, twice. Labelled files present: yes, all three. Section 3 steps 1
to 4, with pure-Python tables: done (`scratch/keep/gate/`). `K_MODEL` chosen and recorded: 3, above.

### Phase 1, the program (2026-10-29 and 30, DONE after three stops; the final state is under ruling 9)

**In short.** The phase stopped three times on its test, each time answered by a ruling, and ends DONE under ruling 9
(the last section of this entry). The final state: `K_MODEL = 3`; the `fcc` suite passes with no frame lost;
71-17/traj01 has its map in the raw-first order. The earlier sections are kept as written at each stop. The program now gates the raw frames before it denoises, in the
order centre, mask, start model, 2D, gates. The lattice gate has its own grade for raw class averages, 0.52, set on
the labelled sets. The model gate's constant, set by ruling 2's rule on the same sets, is 1.5. On the plain simulated
trajectory the model gate at 1.5 removes 9 percent of the frames, where ruling 8 allows 5 percent. 71-17/traj01 went
through the reordered program and has its map: the gates removed nothing, and FSC 0.143 is at 1.12 Å. What follows
first is the entry as written at the first two stops.

The program `trajectory_solve3D` is written and runs: on the real trajectory 71-17/traj01 it
reproduced the previous run's products stage by stage (results below). The test the plan sets for it cannot
pass. The plan asks that the two gates remove at least 90 percent of three windows of twenty pure-noise frames
planted in the simulated trajectory. They removed none of the sixty, and even with the noise planted after the
denoiser they would remove about a quarter. Both gates judge class averages. A class average is a time chunk of
about twenty frames, and a noise window that does not line up with a chunk shares its classes with true frames.
Those classes keep the lattice fringes of their true frames and pass both gates, taking the noise frames with
them. Before the denoiser, the noise frames are moreover replaced by averages of particle frames. Passing would
need the test's bound or the gates' unit (class averages) to change, both fixed by the plan, so the phase stops
here for the maintainer's decision. Everything else of the phase is done and passes: the program, its user
interface, router, parameters, the gate rule with its tester case, the policy note and both test entries'
other checks. Rulings 1 and 3 are applied here for the first time, as written: the new program by that name,
and `autorefine3D_nano` without the atomicity prior with `icm=no` and `lp=1.5`.

**The program as built** (`single_commanders_nano3D.f90`, beside `autorefine3D_nano`).

- Command: `single_exec prg=trajectory_solve3D projfile=<p> stk=<subtr_frames.mrc> smpd=0.358 nthr=24`.
- Stage directories: twelve, `1_trajectory_denoise` to `12_autorefine3D_nano`, each named after the program it
  runs and holding its own project copy. Each stage reads the project file of the stage before and writes its
  own. Stage 5, `refine2D_nano`, reads the import project, as the previous run decided: the ten-chunk
  assignment of stage 3 serves only the diameter.
- Records: `trajectory_solve3D.txt` in the execution directory, written after every stage. One line per stage
  holds its done flag, wall time and numbers: frames and box, the diameters and masks, the classes and the
  fraction of consecutive frames sharing a class, the lattice gate's counts, the model gate's median, spread,
  cut and counts, the shifts, and whether the start model was simulated again. A last line gives the map.
- Resume: a rerun with `mkdir=no` inside the execution directory skips the stages recorded as done. A stage
  that runs again marks every later one as not done.
- Errors: a stage's failure stops the program, which names the stage in its log line
  `>>> TRAJECTORY_SOLVE3D STAGE <n> <program>` and, where the program itself checks, in the error. A stage's
  failure inside the program it calls stops the process with that program's own message; Fortran has no way to
  catch it.
- The model gate (stage 8, inside the commander): `refine3D_nano oritype=cls3D maxits=1 lp=1.0 trs=10 icm=no volrec=no`
  against the start model. `volrec=no` skips a reconstruction from the class averages, which does not change the
  correlations, and the pass's 1.5 GB reprojection models are deleted afterwards. Then
  `model_gate_states(corr, states, K_MODEL)`, then `map_cavgs_selection`, then `prune_particles`. Products: the
  class table `model_gate_classes.tsv` and the class-average/reprojection pairs of the rejected classes and of
  twenty kept ones.
- The rule: `model_gate_states` and `K_MODEL = 3` sit beside `LAT_PEAK_MIN` in `simple_cavg_quality_lattice.f90`,
  and its closed-form case is `test_model_gate` in the lattice tester. The override key is `model_gate_nsd`;
  `simple_parameters_phases.f90` refuses it with any program other than `trajectory_solve3D` and checks the
  program's diameters, `lp` and `thres`. The raw lattice peak in the class table comes from
  `model_cavgs_rejection`, which now writes it as `lattice_peak` on every scored `cls2D` row (a file-table row
  was added for this).

**Decisions where the plan is silent.**

- The diameters are taken at the 0.1 Å `estimate_diam` reports, as the previous run's driver did; that keeps the
  masks identical to the previous run's.
- `element` (default `Pt`) is a key of the start model only.
- `maxits`, `maxits_between` and `nspace` are forwarded to `autorefine3D_nano`. The one-thread test needs them to
  keep its own three rounds of four iterations.
- `autorefine3D_nano` now only defaults to `mkdir=yes` instead of forcing it. The program runs it in its own
  stage directory with `mkdir=no`. Inside `simple_test_exec` no program makes an execution directory, and this
  keeps the layout the same in both executables.
- The atomicity-prior refinement the `fcc` suite exercised (`rec_backend=pcg atm_prior=yes maxits_atm=2`) moved
  to the `wurtzite` suite, which keeps the hand chain. Otherwise the loop would lose its only high-level
  coverage.

**Bugs found and fixed.**

1. A failed run that passed. The parameter layer stopped with a plain `stop` (exit status 0) when an input
   volume or the input stack was missing, so CTest counted the failure as a pass. The first run of the changed
   `fcc` test "passed" in 86 s without making a map. Both checks now use `THROW_HARD` (`error stop 1`) in
   `simple_parameters_phases.f90`. The same pattern remains in `simple_parameters_registry.f90` (a refused file
   format); that file is outside this phase's table and is left open.
2. My own first commander copied numeric keys to the stages as text, so `maxits` arrived as 0. It now uses
   `cmdline%copy_arg`.

**The test as first written** (`single_workflow`, `simple_commanders_test_single.f90`; replaced after ruling 4, below).

- What `fcc` now does: it plants three windows of twenty frames of Gaussian noise, each frame with its own mean
  and spread, drawn from a fixed seed. The windows start at frames 58, 115 and 271 of 300 and are reported in
  the metrics. The suite then runs `trajectory_solve3D` in its stage directory `4_trajectory_solve3D` and
  measures, on the pruned project of the model-gate stage, the noise frames removed (at least 0.90) and the true
  frames lost (at most 0.05). The present floors on the map and the pose checks follow. Pose errors map each
  particle to its frame through `indstk` and leave surviving noise frames out.
- What `wurtzite` does: it keeps the hand chain, with the atomicity-prior refinement.
- Result, `wurtzite` (`scratch/keep/ctest_wurtzite_p1.log`, 902 s with other jobs beside it): passed, with
  the atomicity-prior products present, map correlation 0.953, FSC 0.143 at 1.59 Å, and rotation errors
  2.21° (`autorefine3D_nano`), 2.08° (`refine=cont`) and 2.09° (the polish), within 0.5° of the stage before.
- Result, `fcc` (`scratch/keep/ctest_fcc_p1_2.log`, 467 s): failed on the noise bound alone. Noise removed 0.00
  (bound 0.90), true frames lost 0.00 (bound 0.05), map correlation with the truth 0.985 (floor 0.90), Fourier
  shell correlation (FSC) 0.143 at 3.58 Å (ceiling 5). The pose checks did not run, because the test stops at the
  bound before them.
- What happened inside: refine2D_nano made 15 classes (a fraction 0.85 of consecutive frames share a class), and
  all 15 passed the lattice gate with lattice peaks 0.85 to 1.54. One class was pure noise (15 of 15 frames). It
  passed at 0.85, and its correlation 0.308 sat just above the model gate's cut of 0.306. The other 45 noise
  frames were mixed into classes with true frames (11 of 21, 10 of 30, 9 of 19, ...).
- Why the pure-noise class looks like a particle: the flex denoiser replaces each frame by the pre-image of the
  frames that look like it in a latent model of the whole trajectory. A noise frame comes out as a pooled
  average of particle frames, with their fringes (its record pools a median of 10.2 frames against 7.4 for true
  frames).
- Diagnostic (`scratch/keep/diag/after_denoise/`): the same windows were planted after the denoiser, and the
  program ran stages 2 to 8 on that stack. The lattice gate then rejects the one pure-noise class (16 of 16
  frames), but the other 44 noise frames sit in particle classes (11 of 21, 10 of 30, 6 of 16, ...), and the model
  gate rejects nothing (median 0.79, cut 0.35). That is 0.27 of the noise removed. So the bound cannot be reached
  by class-average gates wherever the noise is planted, as long as the windows straddle the class chunks.

**71-17/traj01 through the program** (`scratch/keep/pipeline/71-17/traj01/`; log
`trajectory_solve3D.log`; driver `scratch/keep/pipeline/run_trajectory.sh`).

- Denoising: the denoiser's table is identical to the previous run's. The program's denoised stack is
  byte-identical to the standalone `trajectory_denoise` of this build, which in phase 0 reproduced the previous
  run's stacks bit for bit; the previous run had deleted this trajectory's stack.
- The diameters: 24.7 / 34.2 Å, giving `mskdiam` 51.3 and `moldiam` 24.7, identical to the previous run's.
- 2D: 649 classes, 45 of them empty, against 649 and 45 in the previous run. The fraction of consecutive frames
  sharing a class is 0.7612 (previous 0.7616).
- The lattice gate: 597 classes accepted, 52 rejected, 12,838 frames kept (previous 598, 51, 12,846; the 2D is
  not bit-repeatable). Run by hand on the program's own 2D project with the same mask, it gives class states
  identical to the program's (same checksum of the `cls2D` states).
- The model gate: median correlation 0.323, robust spread 0.179, cut −0.214. No class was rejected, as the
  labelled sets predicted (phase 0).
- Centering: shift median 9.0 px, largest 43.6 px.
- The 3D: the map was made (`12_autorefine3D_nano/final_results/recvol_state01_iter005.mrc`).
  - The second diameter estimate gives 24.6 / 27.4 Å, so `mskdiam` 41.1 Å (previous run 24.3 / 30.3 Å, 45.5 Å).
    `moldiam` changed by less than 5 percent, so the start model was not simulated again.
  - The 3D took 3,128 s at 24 threads with tests running beside it (previous run 2,103 s), and the whole program
    57 minutes; the largest resident memory was 17.6 GB.
  - The FSC of the final halves at 0.143 is 1.15 Å with the program's mask of 41.1 Å. The previous run reports
    1.36 Å with its 45.5 Å mask, so the two are not like for like.
  - Atoms: 305 strong and 49 weak (previous 364 and 31; the previous run found this count not repeatable from run
    to run, 303 and 82 in its first attempt). The fit correlation is 0.792 (previous 0.775), and the residual
    over noise 1.75 (previous 1.94).

**Files changed in this phase.**

- `src/main/commanders/single/single_commanders_nano3D.f90`: the commander, and `autorefine3D_nano`'s `mkdir`
  default.
- `src/main/ui/single/single_ui_trajectory.f90`, `src/main/exec/single_exec_trajectory.f90`.
- `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90`,
  `simple_parameters_registry.f90` (ruling 4).
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90` (ruling 4).
- `src/main/cavg_quality/simple_cavg_quality_lattice.f90` and its tester.
- `src/main/commanders/simple/simple_commanders_cavgs.f90`.
- `src/main/commanders/test/simple_commanders_test_single.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90`,
  `production/CMakeLists.txt` (comment only), `doc/policies/test_environment_policy.md`.
- `doc/policies/sieving_and_rejection/model_cavgs_rejection.md`, and this note.
- `git diff --summary` shows no mode changes.

**Checks that pass.**

- Release build.
- `single_workflow_wurtzite`.
- The fast gate, 15 of 15 (`scratch/keep/fast_gate_p1_2.log`).
- `unit_single`, `unit_ui`, and `unit_numerics` including the new model-gate case (`scratch/unit/`).
- `check_test_registry.py` and `check_descr.py`.

**After ruling 4 (2026-10-29, late).** Ruling 4 accepted the first stop's findings. The denoiser turns a
structureless frame into a pooled particle, and a class gate cannot remove a frame from a class it keeps. The
ruling replaced the test's junk, kept the model gate as set, and had two defects outside the table fixed. Two
rows were added to the file table for those fixes. Ruling 4 is applied here for the first time.

- The two defects, fixed (fast gate 15 of 15 after both, `scratch/keep/fast_gate_p1_3.log`):
  - `simple_parameters_registry.f90`: a refused file format stopped with exit status 0. It now stops with
    `THROW_HARD`, naming the format, the file and the key; checked by hand with `stk=foo.txt`, which gives exit 1
    and the message.
  - `simple_refine3D_strategy.f90`: `refine3D oritype=cls3D` read past an empty `cls3D` field (a segmentation
    fault in phase 0). It now stops with a message naming the entry count, the class-average count and the
    project, in both the shared-memory and the distributed initialisation (`check_cls3D_field`). Checked on the
    previous run's `make_cavgs` project, which was made before phase 0's fix: exit 1 with "the cls3D field holds 0
    entries for 10 class averages".
- The test's junk, as the ruling sets it. Two windows of sixty frames each (three 2D chunks of twenty) start on a
  chunk boundary, drawn from a fixed seed (they start at frames 41 and 241 of 300). In them, each frame is
  replaced by its clean reprojection, low-passed to 6 Å, normalised and given noise at the trajectory's
  signal-to-noise ratio (SNR) of 0.2, as `stackops snr=` makes the trajectory. The assertions: every junk class
  (at least 80 percent planted frames after the 2D) is removed by the gates, and at most 5 percent of the true
  frames are lost. Reported: the fraction of junk frames removed, the junk frames kept in mixed classes, the
  windows, the seed and where the junk was planted.
- Planted before the denoiser (`scratch/keep/ctest_fcc_p1_3.log`, 346 s; tables in
  `scratch/keep/single_workflow_logs_p1/fcc_junk_before_denoiser/`):
  - The denoiser gave the blob frames lattice: their classes came out with lattice peaks 0.68 to 0.98, above the
    lattice gate's 0.65, against 1.29 and up for the true classes.
  - The 2D mixed blobs with true frames: a fraction of only 0.53 of consecutive frames shared a class.
  - The model gate removed the five blob classes (correlations 0.21 to 0.29 against a median of about 0.85),
    and the lattice gate one more. Both junk classes and every junk frame went, but 28 percent of the true
    frames went with them in the mixed classes (bound 5 percent).
  - By the ruling's own condition (blob classes above 0.65 after denoising), the junk was moved to after the
    denoiser.
- Planted after the denoiser, the form the test now has (`scratch/keep/ctest_fcc_p1_4.log`, 299 s; tables in
  `scratch/keep/single_workflow_logs_p1/fcc_junk_after_denoiser/`):
  - How it runs: the test runs the program's first stage itself, the denoiser into
    `4_trajectory_solve3D/1_trajectory_denoise/`, plants the junk in its output, and records the stage as done
    in `trajectory_solve3D.txt`. The program then resumes from the import with `mkdir=no`, which also exercises
    the resume path.
  - The 2D made five junk classes, each 100 percent planted frames (classes of 15, 33, 17, 30 and 8 frames). The
    lattice gate removed four of them, with lattice peaks 0.46 to 0.64.
  - The fifth, 33 frames, has a lattice peak of 0.70 and a correlation to the start model of 0.26. The model
    gate's cut is −0.25: survivors' median 0.69, robust spread 0.31. The spread is wide because three small
    classes of true frames correlate at only 0.32 to 0.41. So the class passed both gates.
  - Junk classes removed 0.80 (bound 1.0): failed. True frames lost 0.00 (bound 0.05). Junk frames removed 0.58
    in all, and 17 junk frames were kept inside mixed classes. The map correlates 0.983 with the truth, FSC 0.143
    at 3.58 Å; both pass.
- The model gate stays as set (ruling 2's rule, `K_MODEL = 3`), and this is the result, stated plainly. The
  correlation to the start model ranks the lattice gate's survivors well: on the labelled sets the area under the
  curve among them is 0.944, 0.969 and 0.998. But at 98 percent recall, a cut placed relative to the survivors'
  median and spread has no room. It removed 0, 9 and 0 classes on the labelled sets and none on 71-17/traj01.
  In the test it removes nothing either, because true classes of low correlation widen the spread. The program
  delivers what a maintainer needs to look further: the per-class table with correlations and orientations, the
  class-average/reprojection pairs, and the override key `model_gate_nsd`. The remedy for a tail inside accepted
  classes is a frame-level gate after the 3D, outside this run. The per-frame tail measure of section 4 is its
  evidence.
- The test descriptions (`test_environment_policy.md`, the `single_workflow` help text, the CMake comment) state
  the test as it now stands.
- `single_workflow_wurtzite` was run again after the two fixes and passes (`scratch/keep/ctest_wurtzite_p1_2.log`,
  652 s): map correlation 0.954, FSC 0.143 at 1.59 Å, rotation errors 3.15° (`autorefine3D_nano`), 3.17°
  (`refine=cont`) and 3.11° (the polish), each within 0.5° of the stage before. The fast gate (15 of 15),
  `unit_ui`, `unit_single`, the registry and description checks all pass after the last build
  (`scratch/keep/fast_gate_p1_4.log`, `scratch/unit/`).

**What the maintainer needs to decide (second stop).** The test of ruling 4 fails on one pure-junk class that both
fixed gates let through. Its lattice peak is 0.70 against the lattice gate's 0.65, and its model-gate correlation of
0.26 sits above a cut of −0.25. The options seen here:

1. Raise the lattice gate's threshold or tighten the model gate's constant. Both are fixed by rulings, and a
   threshold moved for a test is what ruling 3 excludes.
2. Assert on the junk the gates do catch, for example require the lattice gate to remove the junk classes whose
   lattice peak lies below 0.65 and report the others.
3. Give the model gate an absolute floor on the correlation beside the relative cut. The blob class sits at 0.26
   while every true class of the test is at 0.32 or above, and on 71-17/traj01 the survivors' median is 0.32, so
   such a floor would need setting on the labelled sets as `K_MODEL` was.

The program, its records and resume, the rule as fixed, and the fixes are not affected by the choice.

**Under rulings 5 to 8 (2026-10-30).** The maintainer answered the second stop with four rulings, and each is
applied here for the first time.

- Ruling 5: the test asserts what the gates' policy promises.
- Ruling 6: reject on the raw frames and denoise only what is kept; the gates' constants are set again on raw class
  averages.
- Ruling 7: the order centre, mask, start model, 2D, then the two coupled gates.
- Ruling 8: no planted junk. The `fcc` suite runs the program on the plain trajectory and asserts the map's floors,
  the pose checks and that the gates lose at most 5 percent of the frames.

The phase stopped a third time, on that last bound.

- **The constants on raw class averages** (ruling 6). The maintainer's classes were rebuilt on the raw frames
  (`subtr_frames.mrc`) with the stored assignment and alignment, and phase 0's study repeated
  (`scratch/keep/gate_raw/`: `run_set_raw.sh`, `tables_raw.py`, `lattice_scan.tsv`, `gates*.tsv`, `auc.tsv`).
  - A defect found first: the nano validity gates hard-rejected 55 and 171 non-empty raw class averages on the
    two 71-6 sets, 18 and 43 of them selected by the maintainer. Their "flat image" test measured the local
    variances after the cryo-EM low-pass of 10 Å, and raw graphene-subtracted averages hold almost nothing below
    10 Å. Their spread is the same as that of accepted ones (sd 0.12). The flatness is now measured on the
    unfiltered average (`measure_nano_locvar` in `simple_cavg_quality_feats.f90`), with a tester case: a lattice
    average high-passed at 5 Å must pass. Afterwards only the empty classes are hard-rejected.
  - The lattice peak still ranks the raw averages: the area under the curve (AUC) is 0.990, 0.971 and 0.971,
    pooled 0.978, well above the stop line of 0.85.
  - But 0.65 keeps only 0.984, 0.962 and 0.939 of the selected classes. The rule (the strictest threshold, on a
    0.01 grid, that keeps at least 0.98 on every set) gives 0.52: recall 0.995, 0.995, 0.985, removal of the
    rejected 0.22, 0.20, 0.37; at 0.53, 71-6/traj04 falls to 0.977. Both start models (the diameter on the raw
    labelled classes, and phase 0's chunk diameter) give the same.
  - 0.52 is now `LAT_PEAK_MIN_RAW` beside `LAT_PEAK_MIN` and the program passes it as `thres`. The policy note says
    why there are two.
  - `K_MODEL` by ruling 2's rule on the same tables: behind 0.52 every K in the list keeps the same recall (0.995,
    0.995, 0.985), so the rule's smallest K is 1.5, now `K_MODEL`; on denoised averages it was 3. The model gate
    then removes another 0.02, 0.00 and 0.01 of the rejected classes.
- **The program in the order of ruling 7** (`single_commanders_nano3D.f90`; eleven stages; `trajectory_solve3D.txt`
  records the frames in and out of each):
  1. `trajectory_import` of the raw frames;
  2. `center2D_nano ncls=10`;
  3. `estimate_diam` on the centred chunk averages (`mskdiam = 1.5 x max`, `moldiam = min`, the one estimate);
  4. `simulate_nanoparticle`;
  5. `refine2D_nano mskdiam=<it> center=no`;
  6. `model_cavgs_rejection quality_model=nano thres=0.52` without prune;
  7. the model gate coupled with it, then the prune, each particle's original frame recorded as `frame`, and the
     2D shifts mapped to `ptcl3D`;
  8. `trajectory_extract_substk state=1`;
  9. `trajectory_denoise` (flex);
  10. `trajectory_swap_stack`;
  11. `autorefine3D_nano` with the centred mask and `icm=no`.

  Following `analysis2D_nano`, the refinement starts from the centering's project. The `frame` field is needed
  because `trajectory_swap_stack` renumbers the stack index; the test's pose checks read it, and phase 2's label
  check will.
- **The test** (ruling 8). The planted windows and their code are removed. On the plain simulated trajectory
  (`scratch/keep/ctest_fcc_p1_5.log`, 524 s; tables in `scratch/keep/single_workflow_logs_p1/fcc_rawfirst/`):
  - The centering's shifts are small (median 0.25 px), and the diameter is 17.0 / 17.3 Å, giving `mskdiam` 26.0.
    `refine2D_nano` made 13 non-empty classes (a fraction 0.82 of consecutive frames share a class).
  - The lattice gate at 0.52 kept all 13. The model gate at `K_MODEL = 1.5` removed two of them: median 0.761,
    spread 0.176, cut 0.498, against classes of 11 and 16 frames at correlations 0.427 and 0.460 (lattice peaks
    0.84 and 0.91).
  - Every frame of this trajectory is a true particle, so the gates lost 27 of 300 frames, 0.09 against the bound
    of 0.05: the test fails. The map passes (correlation 0.985, FSC 3.58 Å).
  - Any K above about 1.9 would lose no frame here, and on the labelled sets every K from 1.5 to 4 keeps the same
    recall.
- **71-17/traj01 through the reordered program** (`scratch/keep/pipeline/71-17/traj01/traj01/2_trajectory_solve3D/`,
  22 threads, with the test running beside it; the denoised-first run stays in `1_trajectory_solve3D/`. The log
  file of the same name was overwritten by this run; the earlier run's records, products and the numbers above
  remain):
  - The centering's shifts: median 1.5 px, largest 21.5 px. After `refine2D_nano` with `center=no`: median 2.4,
    largest 36.8. The denoised-first run had 9.0 / 43.6 after its centering, and the previous run a largest of
    82.9.
  - The one diameter estimate is 23.3 / 26.9 Å, so `mskdiam` 40.4 and `moldiam` 23.3.
  - 2D: 576 non-empty classes of 649, a consecutive-same-class fraction of 0.693 (denoised-first 0.761). The raw
    class averages' lattice peaks have a median of 1.22, a 10th percentile of 0.94 and a minimum of 0.53, with 6
    below 0.65; the denoised-first ones had a median of 1.19, 0.90 and 0.51, with 7 below 0.65.
  - The lattice gate at 0.52 rejected only the 73 empty classes. The model gate (median 0.371, spread 0.204, cut
    0.066) rejected none, and so did the coupled gate. All 12,971 frames were denoised (84 s) and refined in 3D.
  - The 3D: map `11_autorefine3D_nano/final_results/recvol_state01_iter005.mrc`.
    - Time: the whole program took 38 minutes and the 3D 2,003 s, at 22 threads with the test beside it; the
      largest resident memory was 17.4 GB.
    - FSC of the final halves: 0.143 at 1.12 Å and 0.5 at 1.85 Å, with the mask of 40.4 Å.
    - Atoms: 374 strong and 25 weak; fit correlation 0.797; residual over noise 1.76.
    - Beside it, the denoised-first run of this phase: 1.15 Å at 41.1, atoms 305 / 49, fit 0.792, residual
      1.75. And the previous run: 1.36 Å at 45.5, atoms 364 / 31, fit 0.775, residual 1.94. The masks differ, so
      the resolutions are not like for like.
- **Other changes and checks.**
  - The policy note gets the two gates of the program in the raw-first order, `LAT_PEAK_MIN_RAW` and why there are
    two grades, `K_MODEL` and the model gate's inertness as a result.
  - The test descriptions (`test_environment_policy.md`, the `single_workflow` help text, the CMake comment) state
    ruling 8's test.
  - The fast gate passes 15 of 15 after the last build (`scratch/keep/fast_gate_p1_6.log`), and the nano tester
    84 of 84 assertions, including the validity case and `K_MODEL`'s closed-form case. `unit_ui`, `unit_single`
    and the registry and description checks pass, and there are no mode changes.
  - `single_workflow_wurtzite` passes after these changes (`scratch/keep/ctest_wurtzite_p1_3.log`, 527 s): map
    correlation 0.980, FSC 0.143 at 2.49 Å, rotation errors 2.19°, 2.05° and 2.01°.

**What the maintainer needs to decide (third stop).** `K_MODEL` set by ruling 2's rule on the raw labelled tables is
1.5, because every K in the list keeps the same recall there and the rule takes the smallest. At 1.5 the model gate
removes two of the 13 classes of the simulated trajectory, whose frames are all true particles: 9 percent against
ruling 8's bound of 5 percent. The rule's constant and the test's bound contradict each other, and choosing another K
is not mine to make. The options seen here:

1. Break the tie on the labelled sets the other way. Where every K keeps the same recall, take the largest (the
   most recall-first, 4), or keep phase 0's 3. Any K above about 1.9 loses no frame of the test.
2. Keep 1.5 and widen the test's bound.
3. Leave the model gate's cut at the override key's discretion and set the default to no cut. The gate removes
   nothing on 71-17/traj01 in either order.

The program, its order, the raw lattice grade, the validity fix and the test's other checks are not affected by
the choice.

**Under ruling 9 (2026-10-30, evening): DONE.** On the raw labelled classes every K keeps the same recall, so ruling
2's rule had nothing to choose with and fell to the strictest value. Ruling 9 breaks such a tie toward recall:
`K_MODEL = 3`, the value phase 0 set on the denoised classes. It is applied here for the first time and recorded in
the constant's comment and in the policy note. The override key stays, and the two gates stay coupled.

- `K_MODEL` is 3 again in `simple_cavg_quality_lattice.f90`. The tester's closed-form case expects 3: cut −0.18956,
  nothing removed. The program's displayed default for `model_gate_nsd` is 3.
- On the raw labelled sets behind 0.52, K = 3 keeps 0.995, 0.995 and 0.985 of the selected classes. It removes
  another 0.02, 0.00 and 0.00 of the rejected ones beyond the lattice gate.
- `single_workflow_fcc` passes (`scratch/keep/ctest_fcc_p1_6.log`, 512 s; tables in
  `scratch/keep/single_workflow_logs_p1/fcc_rawfirst_K3/`).
  - The lattice gate keeps all 13 classes, and the model gate (median 0.761, spread 0.176, cut 0.234) removes none,
    so the frames lost are 0.00 (bound 0.05).
  - The map correlates 0.985 with the truth, FSC 0.143 at 3.58 Å.
  - The rotation errors are 6.37° (`autorefine3D_nano`), 6.08° (`refine=cont`) and 6.04° (the polish), each within
    0.5° of the stage before.
- 71-17/traj01 needs no rerun. Its raw-first run's model gate at 1.5 removed no class (cut 0.066 below every class),
  and at 3 the cut lies lower still, so every stage's decision and the map (1.12 Å, 374 strong atoms) stand.
- Checks after the last build: the fast gate 15 of 15 (`scratch/keep/fast_gate_p1_7.log`), the nano tester,
  `unit_ui` and `unit_single` (`scratch/unit/`), and the registry and description checks all pass. There are no
  mode changes. `single_workflow_wurtzite` passed after the session's last change to its path
  (`scratch/keep/ctest_wurtzite_p1_3.log`).

**Exit items of phase 1.**

- The program, its user interface, router, parameters, the gate rule with its tester case, the policy note, and
  the `single_workflow` change: done, in the raw-first order of rulings 6 to 8.
- The fast gate, both `single_workflow` entries and `unit_single`: all pass.
- One real trajectory, 71-17/traj01, end to end through the program, compared stage by stage with the previous
  run's:
  - In the denoised-first order: the denoised stack byte-identical, the masks identical, the 2D within its
    run-to-run variation, and the lattice gate identical to the hand command.
  - In the raw-first order: the shifts, the lattice-peak distribution, what each gate removed, and the 3D numbers
    beside both. It has its map.
- The files changed over the phase are listed in the sections above. To them add
  `src/main/cavg_quality/simple_cavg_quality_feats.f90` (the nano validity fix, under the `src/main/cavg_quality/`
  row).


### Phase 2, the ten trajectories (2026-10-30, DONE)

**In short.** All ten platinum trajectories went through `trajectory_solve3D` without intervention and have their
maps. The program in its final form (rulings 6 to 9) needed no fix in this phase, so nothing was rerun and no row was
superseded.

- What the gates removed: together they removed 5.6 percent of the frames (7,472 of 132,474). The previous run's
  lattice gate, on denoised class averages, removed 24 percent (31,746). The lattice gate at its raw grade did
  nearly all of the removal: 7,185 frames, between 0 and 140 classes per trajectory. The model gate removed 14
  classes (287 frames) on five trajectories and none on the other five.
- The maps: the resolutions at the program's own masks are the previous run's within 0.04 Å on seven
  trajectories, better on two and worse on one (details below).
- The remaining tail per frame, measured as the plan defines it, is empty on nine trajectories and 4.9 percent on
  one.
- In situ, on the three labelled trajectories, the gates keep 98.5 to 100 percent of the classes that are mostly
  the maintainer's selected frames, and remove 6 to 25 percent of those that are mostly his rejected frames.

**How it was run.** The program was run one trajectory at a time at 24 threads, in the previous run's inventory
order, by `scratch/keep/pipeline/run_queue.sh`:

1. `simple_exec prg=new_project`, then
   `single_exec prg=trajectory_solve3D projfile=<t>.simple stk=/data2/NanoX_trajs/<d>/<t>/subtr_frames.mrc smpd=0.358 nthr=24`
   (`run_trajectory.sh`).
2. Afterwards `post.sh`:
   - computes the Fourier shell correlation (FSC) of the final half maps at the program's mask
     (`11_autorefine3D_nano/fsc_final.log`);
   - applies the keep rule of section 4. The bulky intermediates are dropped: the extracted and the denoised
     stacks (made again in about 80 s by rerunning their stage), the reprojection models, the per-round volumes,
     and the class averages of every 2D iteration but the last. Every project, log, table, montage, triplet stack
     and `final_results/` is kept, 0.4 to 1.2 GB per trajectory, 6.2 GB in all;
   - writes the table row (`table.py`).

71-17/traj01 is the run of phase 1 in the raw-first order. Its model gate then had the constant 1.5 and removed
nothing, so at 3 every decision is the same and it was not rerun. The ten runs took 5.2 hours, 18 to 47 minutes
each.

Disk ran close: during the 3D of the box-200 trajectory (71-17/traj06) free space fell to 16 GB. The extracted stack,
which is no longer needed once denoised, and duplicate copies of the phase 0 and 1 study projects were deleted to
make room. The product of every stage is where the program left it.

**The table** (`scratch/keep/pipeline/table.tsv`, one row per trajectory, with the previous run's values from its
`table.tsv` beside). The columns:

- frames and box, and wall times per stage and in all;
- the centering's and the 2D's shifts;
- the one diameter estimate, the mask and `moldiam`;
- the classes, empty classes, and the fraction of consecutive frames sharing a class;
- what each gate removed (classes and frames), and the model gate's median, spread and cut;
- the frames kept and their fraction;
- the denoiser's choices and the median Kish population (the effective number of frames each frame pools);
- the 3D's low-pass per stage, the FSC at 0.143 and 0.5, the atoms, the fit correlation and the residual over
  noise;
- the per-frame tail: the median and 5th percentile of the 3D's per-frame correlations, the fraction below
  `median - 3 x 1.4826 x MAD`, and the longest run of such frames;
- the space kept.

The main numbers, each followed by the previous run's value in brackets (FSC at the program's mask; the previous
run's masks differ by up to 5 Å, so resolutions are not like for like):

| trajectory | box | lattice gate removed, classes / frames | model gate removed, classes / frames | frames kept [prev] | FSC 0.143, Å [prev] | strong atoms [prev] | largest shift, px [prev] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 71-17/traj01 | 160 | 0 / 0 | 0 / 0 | 12,971 [12,846] | 1.12 [1.36] | 374 [364] | 36.8 [82.9] |
| 71-17/traj02 | 160 | 10 / 129 | 0 / 0 | 12,842 [11,646] | 1.12 [1.12] | 324 [401] | 25.7 [22.1] |
| 71-17/traj03 | 160 | 1 / 2 | 0 / 0 | 12,969 [12,708] | 1.12 [1.12] | 605 [549] | 30.3 [55.7] |
| 71-17/traj04 | 160 | 31 / 445 | 6 / 133 | 12,393 [11,235] | 1.36 [1.36] | 222 [237] | 26.3 [29.1] |
| 71-17/traj05 | 160 | 25 / 302 | 2 / 38 | 12,631 [9,536] | 1.40 [1.36] | 294 [298] | 26.3 [34.9] |
| 71-17/traj06 | 200 | 8 / 87 | 0 / 0 | 12,884 [12,618] | 1.30 [1.30] | 655 [553] | 31.7 [39.7] |
| 71-6/traj01 | 120 | 70 / 1,160 | 0 / 0 | 12,502 [8,186] | 1.13 [1.13] | 224 [209] | 30.2 [35.2] |
| 71-6/traj03 | 120 | 116 / 2,148 | 2 / 35 | 11,479 [4,345] | 1.34 [1.65] | 98 [172] | 25.9 [41.0] |
| 71-6/traj04 | 120 | 140 / 2,811 | 4 / 81 | 10,770 [5,890] | 1.23 [1.26] | 178 [195] | 28.6 [62.6] |
| 71-6/traj02 | 160 | 8 / 101 | 0 / 0 | 13,561 [11,718] | 1.33 [1.10] | 212 [248] | 28.7 [44.6] |

What the table says, in plain words:

- **Far fewer frames are rejected than before.** The raw-first lattice gate at 0.52 removes a small share. On
  71-6/traj01, traj03 and traj04 it removes 8.5, 16 and 21 percent of the frames, where the previous run removed 40,
  68 and 57 percent. Elsewhere it removes 0 to 4.5 percent, against the previous run's 1 to 27 percent. That puts
  2.6 times as many frames into 3D on 71-6/traj03 and 1.8 times as many on 71-6/traj04, and the maps there are no
  worse (1.34 against 1.65 Å, and 1.23 against 1.26 Å).
- **The maps elsewhere are mostly unchanged.** On 71-6/traj02 the program's map is worse (1.33 against 1.10 Å),
  with a mask 1.8 Å larger. On 71-17/traj01 it is better, with a mask 5 Å smaller.
- **The model gate is nearly inert, as phases 0 and 1 found.** Its cut lies below zero or a little above it on
  every trajectory (−0.32 to 0.08), and the survivors' median correlation runs from 0.26 to 0.58.
- **The shifts are tamed.** Centering first and refining with `center=no` keeps the largest shift at 26 to 37 px.
  The previous run's repeated recentering went up to 83 px.
- **The tail per frame is empty.** The 3D's per-frame correlations are broad (median 0.26 to 0.58, robust spread
  comparable to the median), so the plan's cut lies below every frame on nine trajectories. On 71-17/traj03, 4.9
  percent of the frames lie below it, the longest run 23 consecutive frames. A frame-level gate after the 3D, the next
  step named in ruling 4, would need a different measure than this cut.

**In situ on the labelled trajectories** (plan section 3, step 5; `scratch/keep/pipeline/insitu.py`, `insitu.tsv`).
This is a report, not a second setting of a constant. Each class of the workflow's 2D gets the majority label of its
frames (the maintainer's per-frame state in `ptcl2D_field.txt`, through the frame's stack index) and a purity.

- 71-17/traj04: 538 classes of purity 0.8 or more, 84 mostly selected and 454 mostly rejected. The lattice gate
  keeps all 84 and removes 6 percent of the 454; with the model gate, 7 percent. Per frame, 0.8 percent of the
  selected frames and 4.2 percent of the rejected are removed (with the model gate, 0.9 and 5.5 percent).
- 71-6/traj01: 254 pure classes (202 / 52). Recall 0.985, removal 0.25, the model gate adding nothing. Per frame,
  5.0 and 13.9 percent.
- 71-6/traj04: 427 pure classes (101 / 326). Recall 0.990, removal 0.245, or 0.252 with the model gate. Per frame,
  12.5 and 24.8 percent (12.9 and 25.5).

Beside the labelled-set numbers of phase 1 (on the maintainer's own classes rebuilt on the raw frames: recall 0.995,
0.995 and 0.985; removal 0.24, 0.20 and 0.39), the recall holds. The removal is lower on 71-17/traj04, where most of
his rejected frames sit in time chunks that also hold lattice frames. The maintainer selected only 22 percent of that
trajectory's frames, scattered through time.

**Kept for the maintainer's eye** (`scratch/keep/pipeline/<dataset>/<trajectory>/<trajectory>/1_trajectory_solve3D/`;
for 71-17/traj01, `2_trajectory_solve3D/`):

- in `6_model_cavgs_rejection/`: the lattice gate's accepted and rejected class averages
  (`quality_selected_cavgs.mrc`, `quality_rejected_cavgs.mrc`) and its feature table;
- in `7_refine3D_nano/`: the model gate's class table (`model_gate_classes.tsv`) and its class-average/reprojection
  pairs;
- the maps in `11_autorefine3D_nano/final_results/`, with the triplets of class averages against the final map's
  reprojections and the per-frame correlations;
- `trajectory_solve3D.txt`, every stage's numbers;
- `trajectory_solve3D.log` beside the project directory.

**Files changed in this phase.** Only this note. No source changed, so no test or check needed to run again; the
fast gate and both workflow entries last passed at the end of phase 1. The run's scripts are under
`scratch/keep/pipeline/` (`run_queue.sh`, `run_trajectory.sh`, `post.sh`, `table.py`, `insitu.py`); the queue's log
is `scratch/keep/pipeline/queue_phase2.log`.

**Exit item.** Every trajectory mapped: yes, all ten, with exit status 0.

### Phase 3, the report (2026-10-30, DONE)

The report is `doc/implementation_notes/planned/trajectory_solve3D_report.md`. It reads on its own and serves as the
run's review report. It covers:

- what the program does, in its final raw-first order;
- how the gates were set, on denoised and on raw labelled class averages, with the tables' paths;
- what gate 2 removed on the ten trajectories, with the paths of its montages, and the 3D numbers beside the
  previous run's;
- the per-frame tail;
- every fix;
- what changed phase by phase, pointing to `run/review/phase_N.diff`;
- the decisions taken under delegation and the deviations from the plan, with their reasons;
- the evidence, and what is left open.

No source changed in this phase. The files changed are the report and this note. No ruling was added after ruling 9.
