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
