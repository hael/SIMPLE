# The nanoparticle pipeline on every platinum trajectory: report

The unattended run `nano_pipeline_pt` on the Dell, 2026-10-29. The plan, with every detail and number, is
`nano_pipeline_pt.md` beside this report (its section 8, Progress, has one entry per phase). Run files are
under `/home/elmlundho/agent_runs/nano_pipeline_pt/` (below, "the run directory").

## 1. What happened, in short

All ten platinum trajectories in `/data2/NanoX_trajs` went through the present SINGLE pipeline and have a
time-integrated map. No trajectory had to be stopped. The ten are six of dataset `71-17` (12,971 frames each)
and four of `71-6` (13,662 frames each).

The pipeline is the chain of SINGLE programs from tracked frames to a map:
- `trajectory_denoise`, with its new default `pca_mode=flex`, denoises each frame from the frames that look
  like it.
- `refine2D_nano` groups the frames into class averages. Its time-series search keeps neighbouring frames
  together.
- `model_cavgs_rejection quality_model=nano prune=yes` rejects the class averages without a crystal lattice
  and removes their frames.
- `center2D_nano` centres the remaining frames.
- `autorefine3D_nano` refines the 3D map in five stages, with the species-free `detect_atoms` building an
  atomic density model after each stage. It runs without the atomicity prior and with `icm=no`, as ruled.

Two bugs were found, both on the first trajectory, and fixed in the code. The other nine trajectories ran
without intervention through a resumable driver script. The whole run took 4.5 hours of compute at 24
threads.

**What the run did not do.** It did not judge the maps. Ruling 1 reserves that for the maintainer. No map was
compared with the older maps in `/data/NanoX/zenodo`, and the numbers below are what the programs printed, not
quality verdicts.

## 2. What changed in the repository, phase by phase

The diff of each phase is in the run directory as `review/phase_N.diff`.

- **Phase 0** (`review/phase_0.diff`). Both fixes (section 3), the new check in the SINGLE workflow test, and
  the phase 0 Progress entry. Four files changed:
  - `src/main/strategies/parallelization/simple_make_cavgs_strategy.f90` (one line and its comment);
  - `src/main/commanders/single/single_commanders_nano3D.f90` (the transfer table in the copy and delete
    lists);
  - `src/main/commanders/test/simple_commanders_test_single.f90` (one assertion);
  - the plan.
- **Phase 1** (`review/phase_1.diff`) and **Phase 2** (`review/phase_2.diff`). Progress entries only. Eight
  trajectories ran in these phases, none failed, and nothing needed fixing.
- **Phase 3** (`review/phase_3.diff`, written by the driver at the end). This report and the phase 3 Progress
  entry.

No file-table row was added, no test entry was added or registered, no threshold, floor or default was
changed, and no file mode changed.

## 3. The fixes

**Centering after the rejection kept failing on the even/odd halves.** Exposed by `71-17/traj01` at
`center2D_nano`.
- *Symptom.* The run stopped with "particle 4893 has different halves in ptcl2D and ptcl3D". SIMPLE splits the
  frames into two halves, even and odd, so that each half can be reconstructed independently. Both particle
  fields of the project (`ptcl2D` for 2D, `ptcl3D` for 3D) must agree on that split.
- *Cause.* With `tseries=yes`, the class-centre step of the chunk averaging (`make_cavgs`) re-split `ptcl2D` by
  row parity every time it ran. After the rejection has removed frames, rows are no longer frames, so from the
  first removed frame on the new split disagreed with the frames' own halves and with `ptcl3D`.
- *Fix.* Split only when no halves exist yet, the rule the same module already applies in `setup_evenodd`.
  Unpruned projects, including every `analysis2D_nano` run, are split exactly as before.
- *Check.* No tester exists for this module, and the SINGLE workflow test never prunes before centering. The
  check was therefore the rerun on the real trajectory that exposed the bug: it failed before the fix (log in
  `pipeline/71-17/traj01/failed/6_center2D_nano_eo_bug/`) and passed after it, with no frame whose halves differ.
- *Trajectories changed.* Only `71-17/traj01`, which was at that stage; every later trajectory ran with the
  fix.

**`autorefine3D_nano` threw away the transfer table of the atom detection.** Exposed by `71-17/traj01` in 3D.
- *Symptom.* Each round of `detect_atoms` writes `_transfer.tsv`, the map's measured transfer per Fourier
  shell. The detection program's own UI text names it as a product, and the plan lists it among the final
  results. The file was missing from every iteration folder and from `final_results/`.
- *Cause.* The copy and delete lists of `autorefine3D_nano` predate the table, and its final clean-up deleted
  the file.
- *Fix.* The table is added to the copy and delete lists, including the copies made after an atomicity-prior
  stage.
- *Check.* A new assertion in the SINGLE workflow test requires the table in `final_results/`. Both
  `single_workflow_fcc` (1091 s; this suite also runs the atomicity-prior path) and `single_workflow_wurtzite`
  (722 s) pass (`scratch/ctest_single_workflow.log`).
- *Trajectories changed.* `71-17/traj01` was rerun from the 3D stage. Its pre-fix text products are kept in
  `pipeline/71-17/traj01/failed/9_autorefine3D_nano_before_transfer_fix/`. All ten trajectories now have the
  table in `final_results/` and in every iteration folder.

## 4. Evidence against the baseline

The baseline is the maintainer's HEAD `585ce2192` (Release build with tests, gfortran 15.2.1 from
gcc-toolset-15, system FFTW).
- **Fast gate.** 15 of 15 before any change (`scratch/build_clean.log`) and after each fix.
- **Workflow tests.** The two SINGLE workflow tests pass with the new assertion.
- **Fix 1.** The failing real-data stage fails before the fix and passes after it.
- **Repository checks.** `check_test_registry.py` and `check_descr.py` pass.
- **Run-to-run variation.** The 3D refinement does not repeat exactly. `71-17/traj01` ended with 303 strong
  and 82 weak atoms before fix 2, and with 364 strong and 31 weak atoms after it. Fix 2 copies files only, so
  the difference comes from the refinement itself; the second run also used 22 threads instead of 24.
- **Denoiser.** Regenerating a denoised stack gives a bit-identical file.

## 5. The numbers over the ten trajectories

Median, with the range in brackets. The full table is `pipeline/table.tsv` and these figures are in
`pipeline/distributions.txt`, both in the run directory's `scratch/keep/`.

- **Wall time per stage** (seconds, 24 threads):
  - denoise 80 (49–104);
  - the three short stages before the 2D 7 (4–9);
  - `refine2D_nano` 118 (82–171);
  - rejection 1 (0–1);
  - `center2D_nano` 71 (13–126);
  - `autorefine3D_nano` 1370 (456–2254);
  - the whole chain 28 minutes (10–42).

  The 3D is 74 to 90 percent of every trajectory's wall time. The longest 3D runs are on the 160- and 200-pixel boxes with
  many kept frames.
- **Denoising.** Every trajectory got the largest rank offered (32) and `tau` 8. The temporal bandwidth was 8
  frames, except on three `71-6` trajectories where it was 16. z, the denoiser's gain over the parity mean,
  is 91 (28–121): 87–121 on `71-17` and 28–43 on `71-6`.
- **2D.** The fraction of consecutive frames that share a class is 0.75 (0.65–0.86). This is the time-chunk
  behaviour `tseries` is meant to keep; before the fix of 2026-10-28 it was 0.5 percent.
- **Rejection.**
  - Kept fraction of frames: 0.86 (0.32–0.99). On `71-17` it is 0.74–0.99, on `71-6` 0.32–0.86.
  - Classes accepted: 528 (212–604) of 649 or 684.
  - The hard gates reject 31 (6–100) classes. Where checked, every one is a class that `refine2D_nano` left
    empty.
- **Diameters (Å).**
  - First estimate, on unaligned chunk averages: minimum 23.4 (16.2–29.2), maximum 28.5 (23.9–38.5).
  - Second estimate, on centred averages: minimum 20.9 (17.0–27.0), maximum 24.6 (23.0–32.8).
  - Masks: 2D 42.8 (35.8–57.8), 3D 37.0 (34.5–49.2).
  - The rejection's median particle diameter: 24 (22–32).
- **3D.**
  - Resolution at FSC 0.143 between the final half-maps (FSC: Fourier shell correlation): 1.28 Å (1.10–1.65).
  - Atoms: 273 strong (172–553) and 16 weak (3–76).
  - Residual over noise 1.39 (1.11–1.94); `corr_fit` 0.79 (0.75–0.86).
  - Transfer rounds 6 (2–6).

## 6. Where the products are

Everything is under the run directory's `scratch/keep/pipeline/`, one folder per trajectory, with
`README.txt` describing every file and table column. One trajectory's map is
`<dataset>/<trajectory>/<trajectory>/9_autorefine3D_nano/final_results/recvol_state01_iter005.mrc`. Beside it
are the half-maps, the atomic model (`_ATMS.pdb`, `_atoms.csv`) and the density model (`_SIM`, `_FIT`,
`_RES`, `_model.txt`, `_transfer.tsv`). The accepted and rejected class averages for review are in
`5_model_cavgs_rejection/` (`quality_selected_cavgs.mrc`, `quality_rejected_cavgs.mrc`). Every stage keeps its
log, its wall time and its project. The ten trajectories occupy 5.7 GB.

## 7. Decisions taken under delegation

Each is recorded in the plan's Progress entry of the phase named.
- **Run order** (phase 0). Frame count, then box, then name.
- **2D input** (phase 0). `refine2D_nano` reads the import's project, not the chunk-average project, whose
  ten-class assignment serves only the first diameter estimate.
- **Stage layout** (phase 0). Each stage is SIMPLE's own numbered execution directory. `estimate_diam` and the
  denoiser get directories of their own, because they write into the working directory.
- **Montages** (phase 0). They are the rejection program's own `quality_selected_cavgs.mrc` and
  `quality_rejected_cavgs.mrc`, so no `stackops` or `selection` call was needed.
- **Parallelism** (phase 0). One trajectory at a time at 24 threads. The 3D keeps about 20 cores busy and
  dominates the wall time.
- **Keep policy** (phase 0). Keeping every denoised stack would have left about 19.5 GB free, below the plan's
  20 GB. So the stacks of the seven trajectories with box 160 or 200 are not kept, and the driver regenerates
  them bit-identically if a stage must be rerun. The three 120-pixel stacks are kept. Every stage's project is
  kept so that a rerun can start anywhere. The table column `stack_kept` records this per trajectory.
- **Disk** (phase 1). When `/home` fell to 24 GB free during the 200-pixel trajectory, the 5.1 GB of binary
  outputs of the two workflow tests, which had already passed, were deleted and their logs kept.

## 8. Deviations from the plan, and why

Each deviation is the code deciding where the plan's command was refused, as the plan's section 2 provides.
All are recorded in the phase 0 Progress entry.
- `new_project` belongs to `simple_exec`, and `estimate_diam` to `single_exec`.
- `make_cavgs` needs `nparts=1`.
- `center2D_nano` needs `ncls=10`: its UI makes `ncls` required, so the commander's default of ten cannot be
  reached from the command line.
- `autorefine3D_nano` needs `lp`: its UI makes it required, with a documented default of 1.5 Å, which was
  used.

The FSC for the table is computed with `simple_exec prg=fsc` on the final halves, as the phase notes direct.

## 9. Left open for the maintainer

- **The maps themselves.** Their validation is the maintainer's, as ruled.
- **The 2D mask rule on small boxes.** On `71-6/traj03` and `traj04`, 1.5 times the largest diameter of the
  first estimate exceeds the 120-pixel box (43.8 and 43.3 Å against 43.0 Å). The parameter layer then warns and
  uses the box's default mask, 40.8 Å, as designed. The first estimate is made on unaligned chunk averages and
  is larger than the centred one on every trajectory (the 3D mask is 2 to 21 percent smaller than the 2D mask).
  Whether the 2D mask should come from a centred estimate is a choice for the maintainer.
- **Empty classes.** `refine2D_nano` leaves 6 to 100 classes empty after its 15 iterations; the rejection's
  hard gates remove them.
- **Large 2D shifts.** The shifts `refine2D_nano` gives the frames, from class centering repeated over its
  iterations, reach 22 to 83 pixels on some frames, against medians of 3.5 to 8.9. The largest sit on frames
  with little particle signal.
- **Transfer rounds.** The final atom detection reached its cap on transfer rounds on five of the ten
  trajectories (`transfer_rounds_capped = yes` in `_model.txt`).
- **`71-6` behaves differently from `71-17`.** It denoises with a smaller gain and keeps fewer frames
  (0.32–0.86 against 0.74–0.99). Its trajectory `traj03` kept the fewest frames (4,345) and has the coarsest
  half-map resolution (1.65 Å).
- **Repeatability of the 3D.** The 3D refinement is not repeatable run to run (section 4), so two runs of one
  trajectory differ in their atom counts.
- **The UI requirements.** The required `ncls` of `center2D_nano` and `lp` of `autorefine3D_nano` were met,
  not changed. Whether those UI requirements are intended is the maintainer's call.
