# The nanoparticle pipeline on every platinum trajectory: denoise, reject, center, refine

An unattended run on the Dell (`nano_pipeline_pt`), 2026-10-29. It runs the present SINGLE pipeline, from the
graphene-subtracted frames of a tracked particle to its time-integrated map, on every platinum trajectory in
`/data2/NanoX_trajs`, fixes what breaks on the way, and leaves the maps and their products in one place for the
maintainer's own validation. It is a production run with a repair mandate, not a validation study: the run makes
no claim about the maps beyond that the pipeline produced them and what its numbers were.

## 1. What the maintainer fixed (2026-10-29)

- The four stages to exercise, in this order: `trajectory_denoise` with the new default `pca_mode=flex`;
  the automatic class-average rejection (`model_cavgs_rejection quality_model=nano`, the raw lattice gate);
  centering (`center2D_nano`); `autorefine3D_nano` as integrated with the species-free `detect_atoms`
  (each stage refines against the model's `_FIT.mrc`).
- `autorefine3D_nano` runs without the atomicity prior in this first test (no `atm_prior`, no `rec_backend=pcg`)
  and with the ICM regularization off (`icm=no`): the images are already denoised.
- The maps are not judged here. It is sufficient that they are generated; the maintainer validates them
  afterwards. The run collects what that validation will need and the numbers each stage prints.
- Bugs are fixed as they are found, in the phase that finds them, in any file of the chain (section 7), with
  the fast gate run again and the affected trajectories rerun from the fixed stage.
- Every trajectory of `/data2/NanoX_trajs` is run. The data are read only.

Terms. A trajectory is one tracked particle: a stack of frames, one per time step, graphene-subtracted
(`subtr_frames.mrc`). The time-integrated map is the map `autorefine3D_nano` refines from all kept frames of
the trajectory. A class average is the average of the frames `refine2D_nano` groups together; the rejection
works on those. The lattice gate accepts a class average whose strongest Bragg peak stands at least `10^0.65`
times above its ring (`LAT_PEAK_MIN`). `mskdiam` is the diameter of the soft mask in Ångström; `moldiam` the
particle diameter the start volume is simulated with.

## 2. The chain

One directory per trajectory, `<dataset>/<trajectory>/`, under the run's `scratch/keep/pipeline/`. The stages,
each in its own subdirectory, each a `single_exec` or `simple_exec` call with `nthr=24` and `smpd=0.358`
(the header pixel size of these stacks has been wrong before, 0.0358; the key overrides it and the inventory
of phase 0 records what the headers say):

1. **Project.** `single_exec prg=new_project projname=<trajectory>`.
2. **Denoise.** `single_exec prg=trajectory_denoise stk=<data>/subtr_frames.mrc smpd=0.358 nthr=24`, the
   defaults (`pca_mode=flex`, `neigs=32`, `halfwin=0`, `eo_pool=yes`, the inscribed-circle similarity mask).
   Products: `trajectory_denoised.mrcs`, `trajectory_denoise_flex.txt` (per frame: Kish population of the
   weights, the rank, `tau`, `h`), and the log's cross-validation lines (`CV by temporal bandwidth`,
   `CV by latent rank`, fit band, z), which the table keeps.
3. **Import.** `single_exec prg=trajectory_import projfile=<p> stk=trajectory_denoised.mrcs smpd=0.358`.
4. **Diameter, first estimate.** The 2D stages need `mskdiam` before any class average exists, as
   `analysis2D_nano` does it: `simple_exec prg=make_cavgs projfile=<p> ncls=10 tseries=yes refs=chunks.mrc`
   (ten time-chunk averages), then `simple_exec prg=estimate_diam stk=chunks.mrc smpd=0.358 mskdiam=0 lp=3`;
   `mskdiam = 1.5 x max_diam`, `moldiam = min_diam` (the `analysis2D_nano` rule). Recorded.
5. **2D.** `single_exec prg=refine2D_nano projfile=<p> mskdiam=<mskdiam> nthr=24`: its defaults
   (`tseries=yes`, `refine=greedy`, `nptcls_per_cls=20`, `maxits=15`, `center=yes`, `hp=3`, `lp=1`,
   `objfun=cc`, `ml_reg=no`). Products: the class averages of the last iteration and the project's `cls2D` and
   `ptcl2D` fields. The run records the fraction of consecutive frames that share a class (the time-chunk
   behaviour `tseries` is meant to keep; it was 0.5 percent before the fix of 2026-10-28 and should now be
   high); if it is not, that is the first bug to find.
6. **Rejection.** `simple_exec prg=model_cavgs_rejection projfile=<p> quality_model=nano mskdiam=<mskdiam>
   prune=yes nthr=24`. Products: the quality table, the class states, the pruned project (the frames keep
   `indstk`), the median particle diameter over the accepted classes (`diam` on `os_out`), and, for the
   maintainer's review, two montages of the class averages: accepted and rejected (`cavgs_accepted.mrc`,
   `cavgs_rejected.mrc`, by `stackops` or the selection), with the counts and the kept fraction of frames.
7. **Centering.** `single_exec prg=center2D_nano projfile=<p> nthr=24` on the pruned project (it deletes the
   2D solution, keeps the states, makes ten chunk averages, aligns in plane, and carries the shifts to 3D).
   Products: the centered chunk averages and the project.
8. **Diameter, second estimate.** `estimate_diam` on the centered chunk averages as in step 4; the two
   estimates and the rejection's `diam` are recorded side by side. The second estimate sets `mskdiam` and
   `moldiam` for 3D (the first is what the 2D used; a difference above 10 percent is reported, not acted on).
9. **Start volume.** `single_exec prg=simulate_nanoparticle outvol=startvol.mrc smpd=0.358 element=Pt
   moldiam=<moldiam> box=<box> nthr=24` (the element here is the lattice of the reference, as in
   `analysis2D_nano`; detection is species-free).
10. **3D.** `single_exec prg=autorefine3D_nano projfile=<p> vol1=startvol.mrc smpd=0.358 pgrp=c1
    mskdiam=<mskdiam> icm=no nthr=24`: its defaults otherwise (`maxits=5` stages of `maxits_between=10`
    iterations, `lp` of the first stage from the program, `trail_rec=yes`). No `atm_prior`, no `rec_backend`.
    Products: `final_results/` (the map, its halves, `_SIM`, `_FIT`, `_RES`, `_ATMS.pdb`, `_atoms.csv`,
    `_model.txt`, `_transfer.tsv`, `_resid_extrema.csv`), the per-iteration text products
    (`iteration_NN/*_model.txt`, `*_atoms.csv`), the final project.

Where a key above is refused or a stage reads something else than described, the code decides: read the
commander, fix the plan's command in the Progress entry with the reason, go on.

## 3. The data

`/data2/NanoX_trajs/<dataset>/<trajectory>/subtr_frames.mrc` on the Dell: the graphene-subtracted tracked frames
of one platinum particle, 0.358 Å. The same directories may hold older products (`trajectory_denoised.mrcs`
from an earlier denoiser, text fields); the run reads only `subtr_frames.mrc` and writes nothing there.
Phase 0 inventories every trajectory: frames, box, header pixel size, file size, and orders them by frame
count. The old maps of the same particles are in `/data/NanoX/zenodo` (read only); they are not used by the
run.

## 4. What is kept, and the table

Per trajectory, in `scratch/keep/pipeline/<dataset>/<trajectory>/`: the denoised stack, `trajectory_denoise_flex.txt`,
the project after the last stage, the class averages before rejection with the quality table and the two
montages, the centered chunk averages, `startvol.mrc`, `autorefine3D_nano`'s `final_results/` and the
per-iteration text products, every stage's log and wall time. Intermediate volumes and stacks beyond those are
deleted as each trajectory finishes. Phase 0 measures one trajectory's footprint; if the projected total does
not leave 20 GB free on `/home`, the denoised stacks of the largest trajectories are not kept (recorded per
trajectory), the rest is.

`scratch/keep/pipeline/table.tsv`, one row per trajectory: dataset, trajectory, frames, box; denoise: wall time,
chosen rank and bandwidth, median Kish population, frames with a Kish population below 2; 2D: wall time,
classes, consecutive-same-class fraction; rejection: classes accepted and rejected, frames kept and their
fraction, `diam`; the two diameter estimates and the `mskdiam`/`moldiam` used; centering: wall time, median
and maximum shift; 3D: wall time, the per-stage `lp`, the final FSC 0.143 resolution between the halves, the
final `detect_atoms` numbers from `_model.txt` (atoms strong and weak, residual over noise, `corr_fit`,
`transfer_rounds`, `S`, `B_p`, `beta`), and the stage at which the trajectory stopped if it did not finish, with
the error. `README.txt` beside it says what each file is.

## 5. Rules of the run

- Fix as you go. A stage that fails or does something wrong is a bug of this run: find the cause in the code,
  fix it in its owning module (section 7; a file not listed is added with the reason before it is edited),
  give the fix a check in the area's tester where one exists, run the fast gate, and rerun every trajectory
  already past that stage from that stage on. Every fix is a Progress item with the trajectory that exposed
  it. A failure the run cannot fix stops the run with the diagnosis.
- No validation claims. The table reports numbers; the Progress entries and the report say what ran, what
  broke, what was fixed, and what the numbers were. Comparisons to the old maps and judgements of quality are
  the maintainer's.
- One trajectory at a time at 24 threads, smallest first, so that results exist early; two at a time at 12
  threads each only if phase 0 shows a stage that does not scale (the Dell's system FFTW is slow above 16
  threads on some transform sizes; a stage that is unexpectedly slow is timed at 16 threads once).
- The chain is resumable per trajectory and per stage (a stage's directory with its `_FINISHED` marker or
  product is not redone), so a restart after a fix costs only the stages from the fix on.
- Repository rules as always: release 4 keeps no backwards compatibility; `THROW_HARD` arguments never continue
  across lines; no goto; no short-circuit on allocatable elements; dead code met is deleted; new OpenMP code
  is checked at 11 and 24 threads; `git diff --summary` is checked for mode changes; `check_test_registry.py`
  and `check_descr.py` run whenever a test table, router or file header changes. Scripts of the run live in
  `scratch/keep/`, not in the repository. No commits.

## 6. Phases

**Phase 0, one trajectory end to end.** Build the tree on the Dell (Release, tests), the fast gate. Inventory
`/data2/NanoX_trajs`. Run the chain of section 2 on the smallest trajectory, fixing what breaks, and record per
stage the wall time, the products, the disk footprint, and the numbers of section 4; then on the second
smallest without intervention, as the check that the chain runs unattended. Decide the parallelism and the
keep policy from the measurements. Write the per-trajectory driver script (`scratch/keep/run_trajectory.sh`,
resumable per stage) and the table script. Ends DONE when two trajectories have their maps.

**Phase 1, the first half.** The remaining trajectories in order of frame count, the smaller half. Fix as you
go; a fix reruns the affected stages of the trajectories already done, in this phase. The table after each
trajectory. Ends DONE with every trajectory of the half either mapped or stopped with its diagnosis.

**Phase 2, the second half.** The same for the rest.

**Phase 3, the report.** The table complete; the second half of this note's Progress; a report
(`doc/implementation_notes/planned/nano_pipeline_pt_report.md`) that reads on its own: what ran, every fix
with its cause and the trajectory that exposed it, the distribution of each stage's numbers over the
trajectories, where the products are and how the maintainer finds a trajectory's map, and what did not
finish and why. The species status notes are not touched.

## 7. File table

A phase edits only its files; a file the work turns out to need is added here, with its phase and a one-line
reason, before it is edited. Testers (`*_tester.f90`) of a listed Fortran file are covered by its row. The
chain's modules are listed for every phase because the repair mandate reaches all of them.

| File | What changes | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/nano_pipeline_pt.md` | this note: Progress, file-table rows | each |
| `doc/implementation_notes/planned/nano_pipeline_pt_report.md` | the report | 3 |
| `src/main/commanders/single/single_commanders_trajectory.f90` | trajectory_denoise, trajectory_import: fixes | each |
| `src/main/flex/simple_flex_traj_denoise.f90` | the flex denoiser: fixes | each |
| `src/main/flex/simple_flex_cls_kernel.f90` | its kernel: fixes | each |
| `src/main/commanders/single/single_commanders_nano2D.f90` | center2D_nano, refine2D_nano: fixes | each |
| `src/main/commanders/single/single_commanders_nano3D.f90` | autorefine3D_nano, refine3D_nano: fixes | each |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | model_cavgs_rejection: fixes | each |
| `src/main/cavg_quality/` | the nano context and lattice bank: fixes | each |
| `src/main/strategies/search/simple_strategy2D_tseries.f90` | the time-series 2D search: fixes | each |
| `src/main/strategies/search/simple_strategy2D_alloc.f90` | the fresh-start rule of the time series: fixes | each |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | refine2D as the nano programs drive it: fixes | each |
| `src/main/strategies/parallelization/simple_make_cavgs_strategy.f90` | the chunk averages: fixes | each |
| `src/main/nano/simple_nanoparticle.f90` | detect_atoms inside the loop: fixes | each |
| `src/main/nano/simple_nano_density.f90` | the density model: fixes | each |
| `src/main/nano/simple_nanoparticle_utils.f90` | detection utilities: fixes | each |
| `src/main/commanders/simple/simple_commanders_atoms.f90` | detect_atoms commander: fixes | each |
| `src/main/commanders/simple/simple_commanders_sim.f90` | simulate_nanoparticle: fixes | each |
| `src/main/commanders/simple/simple_commanders_imgproc.f90` | estimate_diam: fixes | each |
| `src/main/params/simple_parameters.f90` | a key of the chain: fixes | each |
| `src/main/params/simple_parameters_parse.f90` | the same | each |
| `src/main/params/simple_parameters_phases.f90` | the same | each |
| `src/main/ui/single/` | the UI entries of the chain's programs: fixes | each |
| `src/main/commanders/test/simple_commanders_test_single.f90` | a check for a fix, where the area's test lives here | each |
| `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` | if a fix changes what the policy says | each |

## 8. Progress

(one entry per phase, written by the run)

### Phase 0, one trajectory end to end (2026-10-29, DONE)

**In short.** Both trajectories of the phase now have maps, and the second went through the per-trajectory
driver without intervention. The chain needed two fixes. First, `center2D_nano` stopped on the pruned project
because the chunk averages it makes re-labelled the even/odd halves of the frames by row. Second,
`autorefine3D_nano` dropped the transfer table of the atom detection. Four commands needed a different key or
program from the one the plan gives. Rulings 1 to 3 were applied as written: every trajectory runs,
`autorefine3D_nano` runs without the atomicity prior and with `icm=no`, and the maps are reported with their
numbers but not judged.

**Build.** Release build with tests (`./compile_clean.sh`), using the gfortran 15.2.1 of gcc-toolset-15
(`/opt/rh/gcc-toolset-15/root/usr/bin`), system FFTW and CMake 3.26.5, the same as the maintainer's
`~/src/SIMPLE/build`. Fast gate 15/15 passed after the clean build and after each fix. Every job sources
`scratch/keep/env.sh`.

**Inventory** (`scratch/keep/pipeline/inventory.tsv`). There are ten trajectories in two datasets, and both
datasets have more than one trajectory. `71-17` has six trajectories of 12,971 frames each (`traj01` to `traj05`
at box 160, `traj06` at box 200). `71-6` has four of 13,662 frames each (`traj01`, `traj03` and `traj04` at box
120, `traj02` at box 160). Every MRC header gives 0.358 Å per pixel, the correct value, and `smpd=0.358` is
passed anyway. The stacks are 0.79 to 2.08 GB. The run order is frame count, then box, then name, so the first
trajectory is `71-17/traj01` and the second `71-17/traj02`.

**Commands as run.** Each stage runs in its own numbered directory under the project directory
`<dataset>/<trajectory>/<trajectory>/`. SIMPLE makes these directories itself (`N_<program>`, the next free
integer prefix), and every stage reads the project of the stage before it. The denoiser runs outside the
project, in `01_denoise/`. Where the plan's command was refused, the code decided, as section 2 provides:

1. `simple_exec prg=new_project projname=<trajectory>`. The program belongs to `simple_exec`, not `single_exec`.
2. `single_exec prg=trajectory_denoise stk=/data2/NanoX_trajs/<d>/<t>/subtr_frames.mrc smpd=0.358 nthr=24`, as planned.
3. `single_exec prg=trajectory_import projfile=<t>.simple stk=<...>/01_denoise/trajectory_denoised.mrcs smpd=0.358`.
4. `simple_exec prg=make_cavgs projfile=1_trajectory_import/<t>.simple ncls=10 tseries=yes refs=chunks.mrc nparts=1 nthr=24`.
   Without `nparts` the program refuses to run ("Missing key on command line: nparts"); `center2D_nano` passes
   `nparts=1` the same way. Then `single_exec prg=estimate_diam stk=2_make_cavgs/chunks.mrc smpd=0.358 mskdiam=0 lp=3 nthr=24`
   (a `single_exec` program). Without a project it writes into the working directory, so it gets its own
   directory, `3_estimate_diam/`. `mskdiam = 1.5 x max diameter`, rounded to 0.1 Å.
5. `single_exec prg=refine2D_nano projfile=1_trajectory_import/<t>.simple mskdiam=<mskdiam> nthr=24`. It reads
   the import's project, not the one from `make_cavgs`, because the ten-chunk assignment there serves only the
   diameter. `refine2D_nano` replaces the class assignment with its own 20-frame chunks in any case.
6. `simple_exec prg=model_cavgs_rejection projfile=4_refine2D_nano/<t>.simple quality_model=nano mskdiam=<mskdiam> prune=yes nthr=24`,
   as planned. The program writes the two montages itself, `quality_selected_cavgs.mrc` (accepted) and
   `quality_rejected_cavgs.mrc` (rejected), so no `stackops` or `selection` call is needed.
7. `single_exec prg=center2D_nano projfile=5_model_cavgs_rejection/<t>.simple ncls=10 nthr=24`. The program's
   user interface (UI) marks `ncls` as required, so the commander's own default of ten classes cannot be reached
   from the command line.
8. `estimate_diam` as in step 4, on `6_center2D_nano/cavgs_iter<last>.mrc`, in `7_estimate_diam/`.
9. `single_exec prg=simulate_nanoparticle outvol=startvol.mrc smpd=0.358 element=Pt moldiam=<min diameter> box=<box> nthr=24`, as planned.
10. `single_exec prg=autorefine3D_nano projfile=6_center2D_nano/<t>.simple vol1=8_simulate_nanoparticle/startvol.mrc smpd=0.358 pgrp=c1 mskdiam=<mskdiam 3D> lp=1.5 icm=no nthr=24`.
    The UI makes `lp` required, with a documented default of 1.5 Å (the SINGLE workflow test uses the same
    value). The first stage matches at 1.507 Å and the program's later stages at 1.005 Å.
11. For the table, `simple_exec prg=fsc vol1=<final _even> vol2=<final _odd> smpd=0.358 mskdiam=<mskdiam 3D> nthr=24`
    in the 3D directory (`fsc_final.log`).

**Fix 1: chunk averages after a prune keep the even/odd halves** (`simple_make_cavgs_strategy.f90`). Exposed
by `71-17/traj01` at centering. `center2D_nano` stopped with "inconsistent even/odd partition in the project
... particle 4893 has different halves in ptcl2D and ptcl3D". The cause: with `tseries=yes` the class-centre
step of `make_cavgs` re-partitioned `ptcl2D` by row parity every time it ran. After the rejection has pruned
frames, rows are no longer frames. So from the first removed frame (row 4893) on, the new halves disagreed with
the frames' own halves and with `ptcl3D`, which kept them. The halves check of `refine2D` then stopped the run.
The fix partitions only when the field has no halves yet, which is the rule `setup_evenodd` in the same module
already applies. A fresh project (the chunk averages of step 4, and every `analysis2D_nano` run) gets the same
row-parity halves as before. No tester exists for the `make_cavgs` strategy, and the SINGLE workflow test never
prunes before centering. The check is therefore the rerun on the trajectory that exposed the bug: before the
fix it stopped (log kept in `71-17/traj01/failed/6_center2D_nano_eo_bug/`); after it, centering finished with no
frame whose halves differ between the two fields. The fix came before any later stage had run, so no rerun of
other trajectories was needed.

**Fix 2: `autorefine3D_nano` keeps the transfer table** (`single_commanders_nano3D.f90`). Exposed by
`71-17/traj01` in 3D. Each round of `detect_atoms` writes `_transfer.tsv` (the map's transfer per Fourier
shell), a product named in its UI text and listed by this plan for `final_results/`. The copy and delete lists
of `autorefine3D_nano` predate it: the table was never copied into `iteration_NN/` or `final_results/`, and the
final clean-up (`rm recvol*`) deleted it. The fix adds it to the copy and delete lists of each round and to the
copies of the atomicity-prior stage's `_base` products. The check is a new assertion in the SINGLE workflow
test (`simple_commanders_test_single.f90`) that `final_results/` holds the transfer table; both
`single_workflow_fcc` (1091 s, which also runs the atomicity-prior path) and `single_workflow_wurtzite` (722 s)
pass (`scratch/ctest_single_workflow.log`). `71-17/traj01` was rerun from `autorefine3D_nano` (its pre-fix text
products are kept in `71-17/traj01/failed/9_autorefine3D_nano_before_transfer_fix/`). `71-17/traj02` ran with
the fix from the start. Both have `_transfer.tsv` in `final_results/` and in every `iteration_NN/`.

**Observations, not acted on.**
- The time-series behaviour of `refine2D_nano` holds. The fraction of consecutive frames that share a class is
  0.762 and 0.738, against 0.5 percent in the labelled runs before the fix of 2026-10-28.
- The shifts `refine2D_nano` gives the frames (`center=yes`, class centering repeated over 15 iterations) reach
  83 pixels on `traj01`, against a median of 8.9. For example, frames 516 to 555 sit near 60 pixels with
  correlations around 0.05, which suggests frames carrying little particle signal. `center2D_nano` adds only
  small increments to these shifts.
- The 3D refinement does not repeat exactly between runs. The pre-fix and post-fix runs of `traj01` (24 threads,
  then 22 threads beside the two workflow tests) ended with 303 strong and 82 weak atoms, and with 364 strong
  and 31 weak atoms. The fix copies files only, so the difference comes from the refinement itself.
- The rejection's hard gates remove the classes that `refine2D_nano` left empty: 45 of 649 classes on `traj01`
  and 63 on `traj02`. The lattice gate rejects 6 and 65 more.
- On `traj01` the 3D mask is 11 percent smaller than the 2D one (45.5 against 51.3 Å). This exceeds the
  10 percent the plan asks to report, and it was not acted on.

**Wall time, CPU and memory per stage** (`71-17/traj02`, unattended at 24 threads, nothing else running;
CPU use and peak memory from `/usr/bin/time -v` at the end of each stage log): denoise 79 s, 2225 % CPU, 5.1 GB;
import 0 s; `make_cavgs` 7 s; `estimate_diam` 0 s; `refine2D_nano` 118 s, 1767 %, 4.2 GB; rejection 1 s;
`center2D_nano` 83 s, 434 %, 2.4 GB; second `estimate_diam` 0 s; `simulate_nanoparticle` 0 s;
`autorefine3D_nano` 1498 s, 1999 %, 11.2 GB. The whole chain takes about 30 minutes. `71-17/traj01` took 80,
116, 81 and 1930 s for denoise, 2D, centering and 3D (the 3D at 24 threads, the first run), and 2103 s for the
3D rerun at 22 threads beside the tests. Regenerating a denoised stack gives a bit-identical file in 80 s.

**Parallelism.** One trajectory at a time at 24 threads, as section 5 sets by default. The 3D takes about 85
percent of a trajectory's wall time and keeps about 20 of the 24 cores busy. The one stage that scales poorly,
`center2D_nano` at about 4 cores, lasts 83 s. Two trajectories at 12 threads would save at most a few
minutes per trajectory, so no 16-thread timing was needed.

**Disk and keep policy.** Before the keep policy, one trajectory occupies 8.1 GB at its peak (2.4 GB of it the
reprojection models of the 3D). After it, it occupies 1.7 GB, of which 1.33 GB is the denoised stack.
Keeping every stack would leave about 19.5 GB free on `/home` at the end of the run, and less while a 3D is
running. That is below the 20 GB the plan requires. So, as section 4 provides, the denoised stacks of the
larger trajectories (box 160 and 200, seven of the ten) are not kept, and the three 120-pixel stacks of `71-6`
are. This is recorded per trajectory in the table column `stack_kept`. The projected free space at the end is
then about 30 GB. A dropped stack is regenerated bit-identically by the driver whenever a stage of that
trajectory has to run again; this path was exercised once on `traj02`. Each kept trajectory takes about
0.4 GB (`traj01` 396 MB, `traj02` 374 MB). Besides the plan's list, the policy deletes `cavg_state.bin` (the
200 MB of class sums needed only to continue a 2D run), the even/odd halves and per-iteration images of the 2D
class averages, and every volume in the 3D directory outside `final_results/`. It keeps every stage's project
(about 7 MB each), so a rerun after a fix can start at any stage.

**Driver and table.** `scratch/keep/run_trajectory.sh <dataset> <trajectory>` has one function per stage. It
skips a stage whose product and normal-stop line exist, and removes and reruns one that is unfinished. It writes
each stage's log, `walltime.txt` and the `time -v` summary into the stage's directory, stops the trajectory with
`status.txt` (`STOPPED <stage>: <error>`) on a failure, applies the keep policy, and appends the trajectory's row
with `scratch/keep/table.py`. `table.py` replaces an earlier row of the same trajectory and moves the old row to
`table_superseded.tsv`. `scratch/keep/pipeline/README.txt` says what each file and column is.

**The first two rows** (`scratch/keep/pipeline/table.tsv`):

| | `71-17/traj01` | `71-17/traj02` |
| --- | --- | --- |
| frames, box | 12,971, 160 | 12,971, 160 |
| denoise: wall, rank / tau / h, z | 80 s, 32 / 8 / 8 frames, 120.7 | 79 s, 32 / 8 / 8 frames, 96.0 |
| median Kish population, frames below 2 | 8.72, 11 | 10.23, 0 |
| 2D: wall, classes, consecutive same class | 116 s, 649, 0.762 | 118 s, 649, 0.738 |
| rejection: classes accepted / rejected (empty) | 598 / 51 (45) | 521 / 128 (63) |
| frames kept, fraction, `diam` | 12,846, 0.990, 26.5 Å | 11,646, 0.898, 26.0 Å |
| diameters, first: min / max, `mskdiam` 2D | 24.7 / 34.2, 51.3 Å | 24.3 / 27.5, 41.2 Å |
| diameters, second: min / max, `mskdiam` 3D, `moldiam` | 24.3 / 30.3, 45.5, 24.3 Å | 23.7 / 25.6, 38.4, 23.7 Å |
| centering: wall, median / max shift | 81 s, 8.90 / 82.9 px | 83 s, 3.46 / 22.1 px |
| 3D: wall, lp per stage | 2103 s (rerun), 1.507 / 1.005 x4 Å | 1498 s, 1.507 / 1.005 x4 Å |
| FSC 0.143 between the final halves | 1.36 Å | 1.12 Å |
| atoms strong / weak | 364 / 31 | 401 / 3 |
| residual over noise, `corr_fit` | 1.94, 0.775 | 1.40, 0.852 |
| transfer rounds, `S`, `B_p`, `beta` | 4, 150.0, 31.1, -0.054 | 6, 494.5, 62.3, -0.039 |
| kept on disk, stack kept | 396 MB, no | 374 MB, no |

**Exit items.** The tree is built and the fast gate passes: done, 15/15 after each fix
(`scratch/build_clean.log`). The inventory is written: done. The smallest trajectory ran with fixes, and its
wall times, products, footprint and numbers are recorded: done, above. The second smallest ran without
intervention: done, `71-17/traj02/driver.log` and status `DONE`. The parallelism and keep policy are decided:
done, above. The driver and table scripts are written: done. Two trajectories have their maps: done, in
`<dataset>/<trajectory>/<trajectory>/9_autorefine3D_nano/final_results/`. Every path here is under
`/home/elmlundho/agent_runs/nano_pipeline_pt/scratch/keep/pipeline/` unless written otherwise.

### Phase 1, the first half (2026-10-29, DONE)

**In short.** All four trajectories of this half have maps, and none needed a fix. The half is the four
smallest of the eight that remained after phase 0, in inventory order (frame count, then box, then name):
`71-17/traj03`, `traj04`, `traj05` and `traj06`. Each ran unattended through `scratch/keep/run_trajectory.sh`,
one after the other at 24 threads, between 13:20 and 15:33 (`scratch/keep/queue_phase1.log`). No stage failed.
Every log is free of errors, and the only warning is the expected one from the stages that start without a
mask (`estimate_diam` and the inner 2D run of `center2D_nano`, both with `mskdiam=0`). No repository file
changed in this phase apart from this note. The commands are those of the phase 0 entry. Rulings 1 to 3 hold
as there: no atomicity prior, `icm=no`, and no judgement of the maps.

**What the numbers say, per stage, over the four** (rows in `scratch/keep/pipeline/table.tsv`):
- **Denoise.** 78 to 104 s; the 200-pixel `traj06` was the slowest. On every trajectory the cross-validation
  chose the largest rank offered (32), `tau` 8 and a temporal bandwidth of 8 frames, with z from 87 to 102. The
  median Kish population (the effective number of frames each denoised frame averages) is 10.3 to 11.9, and no
  frame falls below 2.
- **2D (`refine2D_nano`).** 117 to 171 s, 649 classes each. The fraction of consecutive frames that share a
  class is 0.738 to 0.848.
- **Rejection.** 1 s or less. Between 0.735 (`traj05`, 464 of 649 classes accepted) and 0.980 (`traj03`) of the
  frames are kept. The hard gates reject 10 to 100 classes. On `traj03`, the 100 hard-gate rejections were
  checked one by one, and every one is a class that `refine2D_nano` left empty. The median particle diameter
  over the accepted classes is 23.5 to 32.0 Å.
- **Diameters.** On every trajectory the 3D mask is smaller than the 2D one, by 5 percent (`traj03`) to
  15 percent (`traj06`, 57.8 to 49.2 Å). `traj04` (−10.9 percent) and `traj06` exceed the 10 percent the plan
  asks to report; this was not acted on.
- **Centering.** 58 to 126 s. The median in-plane shift the 3D receives is 4.4 to 8.1 pixels and the maximum
  29 to 56 pixels; as in phase 0, these come mostly from the class centering of `refine2D_nano`.
- **3D (`autorefine3D_nano`).** 1074 to 2254 s. The matching low-pass limit is 1.507 Å in the first stage and
  1.005 Å after; on the 200-pixel box it is 1.492 and 0.994 Å, the nearest Fourier shells. The resolution at
  FSC 0.143 between the final half-maps is 1.12 to 1.36 Å. The final atom detection found 237 to 553 strong
  and 4 to 76 weak atoms, with a residual over noise of 1.11 to 1.81, `corr_fit` 0.75 to 0.86, and 2 to 6
  transfer rounds.

| | `traj03` | `traj04` | `traj05` | `traj06` |
| --- | --- | --- | --- | --- |
| box | 160 | 160 | 160 | 200 |
| wall: denoise / 2D / centering / 3D (s) | 78 / 118 / 60 / 2254 | 84 / 123 / 88 / 1248 | 85 / 117 / 58 / 1074 | 104 / 171 / 126 / 2084 |
| z, median Kish population | 86.5, 11.87 | 98.5, 10.48 | 102.2, 10.63 | 97.3, 10.30 |
| consecutive same class | 0.738 | 0.839 | 0.848 | 0.765 |
| classes accepted / rejected (hard gate) | 536 / 113 (100) | 548 / 101 (10) | 464 / 185 (13) | 604 / 45 (27) |
| frames kept, fraction, `diam` | 12,708, 0.980, 31.0 Å | 11,235, 0.866, 23.5 Å | 9,536, 0.735, 24.5 Å | 12,618, 0.973, 32.0 Å |
| `mskdiam` 2D / 3D, `moldiam` (Å) | 49.2 / 46.7, 27.0 | 41.1 / 36.6, 20.7 | 39.8 / 37.3, 21.1 | 57.8 / 49.2, 25.0 |
| shift median / max (px) | 6.87 / 55.7 | 4.37 / 29.1 | 4.91 / 34.9 | 8.07 / 39.7 |
| FSC 0.143 | 1.12 Å | 1.36 Å | 1.36 Å | 1.30 Å |
| atoms strong / weak | 549 / 51 | 237 / 76 | 298 / 10 | 553 / 4 |
| residual over noise, `corr_fit` | 1.81, 0.762 | 1.11, 0.856 | 1.16, 0.832 | 1.38, 0.753 |
| transfer rounds, `S`, `B_p`, `beta` | 6, 80.0, 16.4, −0.085 | 2, 1333, 98.8, −0.146 | 6, 102.5, 32.7, −0.134 | 6, 54.6, 15.8, −0.030 |
| kept on disk | 377 MB | 374 MB | 369 MB | 642 MB |

**Disk.** While the 3D of the 200-pixel `traj06` ran, the free space on `/home` fell to 24 GB, below the
30 GB at which the phase notes ask for a stricter keep policy. The finished trajectories' denoised stacks were
already gone (phase 0's policy), so the space came from elsewhere. 5.1 GB were the binary outputs of the two
SINGLE workflow tests that passed in phase 0 (`repo/build/test_runs/single_workflow_*`); these were deleted and
their logs kept. Free space then stayed between 26 and 29 GB until `traj06` finished, and was 37 GB at the end
of the phase. A 200-pixel trajectory occupies about 9 GB while its 3D runs (each reprojection model is 0.9 GB)
and 642 MB once it is finished.

**Exit items.** Every trajectory of the half is mapped: done, four of four with status `DONE` and a map in
`<dataset>/<trajectory>/<trajectory>/9_autorefine3D_nano/final_results/`. The table was updated after each
trajectory: done. Fixes would have rerun the affected stages: none were needed. Paths are under
`/home/elmlundho/agent_runs/nano_pipeline_pt/scratch/keep/pipeline/`.

### Phase 2, the second half (2026-10-29, DONE)

**In short.** The last four trajectories, all of dataset `71-6`, have maps, so every platinum trajectory in
`/data2/NanoX_trajs` is now mapped. As in phase 1, no stage failed and no fix was needed. No trajectory of
phase 1 had to be rerun, and no repository file changed in this phase apart from this note. The trajectories
ran unattended through `scratch/keep/run_trajectory.sh` in inventory order, between 15:34 and 16:45
(`scratch/keep/queue_phase2.log`): `traj01`, `traj03` and `traj04` at box 120, then `traj02` at box 160.
Commands and rulings are as in the phase 0 entry.

**One new warning, by design.** On `71-6/traj03` and `traj04` the 2D mask from the plan's rule (1.5 times
the largest diameter of the first estimate) came out at 43.8 and 43.3 Å. That is wider than the 120-pixel box
(43.0 Å). The parameter layer (`simple_parameters_phases.f90`) handles this case by design: it warns "Mask
diameter too large, falling back on default value" and uses the box's default mask, (120 − 6) × 0.358 =
40.8 Å. This happened three times in `refine2D_nano` and once in the rejection. The code did what it is meant
to do, so this is not a fix case, and the plan's rule was not changed. The first estimate is made on unaligned
time-chunk averages and comes out larger here: its largest diameters were 29.2 and 28.9 Å, against 23.0 and
24.4 Å for the centred estimate that sets the 3D mask (34.5 and 36.6 Å, well inside the box).

**What the numbers say over the four** (rows in `scratch/keep/pipeline/table.tsv`):
- **Denoise.** The denoiser behaves differently on this dataset from `71-17`. On three of the four
  trajectories the cross-validation picks a temporal bandwidth of 16 frames instead of 8, still at rank 32 and
  `tau` 8. The gain over the parity mean (z) is 28 to 43, against 87 to 120 on `71-17`. The median Kish
  population is 19 to 23 on the three bandwidth-16 trajectories and 10.7 on `traj03`. Wall time is 49 to 53 s
  at box 120 and 84 s at box 160.
- **2D.** 82 to 130 s, 684 classes each. The fraction of consecutive frames that share a class is 0.650 to
  0.863.
- **Rejection.** It keeps fewer frames on this dataset than on `71-17`: 0.318 (`traj03`, 212 of 684 classes
  accepted), 0.431 (`traj04`), 0.599 (`traj01`) and 0.858 (`traj02`). The hard gates, which on the checked
  trajectories reject only empty classes, account for 6 to 46 classes; the lattice gate rejects the rest. The
  median particle diameter over the accepted classes is 22.0 to 23.5 Å.
- **Diameters.** The 3D mask is 2 to 21 percent smaller than the 2D mask actually requested. On `traj03`,
  `traj04` and `traj02` the difference exceeds the 10 percent the plan asks to report; this was not acted on.
- **Centering.** 13 to 89 s. The median shift is 3.7 to 8.8 pixels and the maximum 35 to 63 pixels.
- **3D.** 456 to 1595 s. At box 120 the matching low-pass limit is 1.481 Å, then 0.999 Å (the nearest Fourier
  shells). The resolution at FSC 0.143 between the final half-maps is 1.10 to 1.65 Å; the largest value is on
  `traj03`, the trajectory with the fewest kept frames (4,345). The final atom detection found 172 to 248
  strong and 6 to 21 weak atoms.

| | `traj01` | `traj03` | `traj04` | `traj02` |
| --- | --- | --- | --- | --- |
| box | 120 | 120 | 120 | 160 |
| wall: denoise / 2D / centering / 3D (s) | 49 / 82 / 29 / 820 | 53 / 93 / 13 / 456 | 53 / 90 / 18 / 575 | 84 / 130 / 89 / 1595 |
| denoise h, z, median Kish population (frames below 2) | 16, 34.2, 20.7 (6) | 8, 42.6, 10.7 (0) | 16, 28.1, 22.9 (0) | 16, 40.6, 19.1 (8) |
| consecutive same class | 0.650 | 0.863 | 0.678 | 0.716 |
| classes accepted / rejected (hard gate) | 375 / 309 (34) | 212 / 472 (6) | 290 / 394 (19) | 541 / 143 (46) |
| frames kept, fraction, `diam` | 8,186, 0.599, 22.0 Å | 4,345, 0.318, 22.0 Å | 5,890, 0.431, 23.0 Å | 11,718, 0.858, 23.5 Å |
| `mskdiam` 2D requested / 3D, `moldiam` (Å) | 35.8 / 35.1, 18.8 | 43.8 (box default 40.8 used) / 34.5, 17.7 | 43.3 (40.8 used) / 36.6, 17.0 | 42.2 / 35.1, 20.3 |
| shift median / max (px) | 3.73 / 35.2 | 4.80 / 41.0 | 8.77 / 62.6 | 4.74 / 44.6 |
| FSC 0.143 | 1.13 Å | 1.65 Å | 1.26 Å | 1.10 Å |
| atoms strong / weak | 209 / 21 | 172 / 8 | 195 / 6 | 248 / 21 |
| residual over noise, `corr_fit` | 1.83, 0.812 | 1.32, 0.765 | 1.21, 0.745 | 1.65, 0.835 |
| transfer rounds, `S`, `B_p`, `beta` | 6, 27.9, 1.15, 0.030 | 6, 39.2, 8.85, −0.118 | 4, 61.4, 20.7, −0.060 | 6, 51.0, 2.33, 0.043 |
| kept on disk, denoised stack kept | 964 MB, yes | 950 MB, yes | 954 MB, yes | 392 MB, no |

**Disk.** Free space on `/home` stayed between 27 and 36 GB, with 34 GB free at the end. The three 120-pixel
trajectories keep their denoised stacks, as phase 0 decided; the 160-pixel `traj02` does not.

**Exit items.** Every trajectory of the second half is mapped: done, four of four with status `DONE`. Each of
the ten trajectories of the run has `final_results/recvol_state01_iter005.mrc` and its `_transfer.tsv`,
checked directory by directory. The table was updated after each trajectory: done, ten rows. No fix was made
in this phase, so no earlier trajectory needed a rerun. Paths are under
`/home/elmlundho/agent_runs/nano_pipeline_pt/scratch/keep/pipeline/`.

### Phase 3, the report (2026-10-29, DONE)

**In short.** The table is complete and agrees with the directories, and the report is written:
`nano_pipeline_pt_report.md` beside this note. All ten trajectories have status `DONE`, a row in the table and
every product the plan lists. No trajectory stopped, so no row carries an error. The run never judged the
maps, and the report says so.

**The check.** `scratch/keep/distributions.py` reads the table, the inventory and each trajectory's directory.
Its output is `scratch/keep/pipeline/distributions.txt`. It confirms, for each of the ten trajectories:
- one row in `table.tsv` and status `DONE`;
- every final product in `final_results/`: the map, its halves, `_SIM`, `_FIT`, `_RES`, `_ATMS.pdb`,
  `_atoms.csv`, `_model.txt`, `_transfer.tsv` and `_resid_extrema.csv`;
- five per-iteration `_model.txt` files;
- both montages of the rejection.

The same file holds the median and range of every number in the report, plus the transfer-round flags of each
final model. Five of the ten final detections reached the cap on transfer rounds; the report leaves this open
for the maintainer. `scratch/keep/pipeline/README.txt` now also says how to find one trajectory's map, which
denoised stacks are kept, and what the scripts do. The ten trajectories occupy 5.7 GB.

**Files changed in this phase.** `doc/implementation_notes/planned/nano_pipeline_pt_report.md` (new) and this
note. No source file changed, so no build or test was needed. Every job of the run has ended.

**Exit items.** The table is complete and checked against the directories: done (`distributions.txt`). The
second half of Progress is written: done (the entries of phases 1 to 3). The report covers what ran, every
fix with its cause, its check and the trajectory that exposed it, the distribution of each stage's numbers,
where the products are, and what did not finish: done; nothing failed to finish. The species status notes
were not touched.
