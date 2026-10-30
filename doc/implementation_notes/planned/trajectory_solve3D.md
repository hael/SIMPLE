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

## 8. Progress

(one entry per phase, written by the run)
