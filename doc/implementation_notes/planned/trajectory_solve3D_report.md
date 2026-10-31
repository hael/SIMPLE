# `trajectory_solve3D`: report of the run (2026-10-29 to 30)

The run made the SINGLE nanoparticle pipeline one program and gave it a second rejection of class averages. It set
that gate's constant on the maintainer's three labelled trajectories, and ran all ten platinum trajectories of
`/data2/NanoX_trajs` through the program, beside the previous run (`nano_pipeline_pt`), which ran the same stages
by hand. The plan is `trajectory_solve3D.md` in this directory; its section 8 (Progress) holds the detail behind
every number here. Paths of products are under `/home/elmlundho/agent_runs/trajectory_solve3D/` (written `run/`
below). The maps are tabled, not judged: the maintainer validates them.

## What the program does

`single_exec prg=trajectory_solve3D projfile=<p> stk=<subtr_frames.mrc> smpd=0.358 nthr=24` takes the
graphene-subtracted frames of one tracked particle to its time-integrated map. Each stage runs in its own numbered
directory, on its own copy of the project. Its final order, set by the maintainer's rulings 6 and 7 during phase 1,
differs from the plan's section 2. Nothing is rejected on evidence the denoiser has touched, because the flex
denoiser turns a frame without a particle into a pooled particle image (measured twice in phase 1).

1. Import of the raw frames.
2. `center2D_nano` on ten time-chunk averages: in-plane alignment, the shifts carried to the 3D field.
3. `estimate_diam` on the centred averages, the one diameter estimate of the run: `mskdiam = 1.5 x` the largest
   diameter, `moldiam` the smallest.
4. `simulate_nanoparticle`: the start model, a face-centred cubic (FCC) platinum particle of that diameter.
5. `refine2D_nano` on the raw frames with that mask and `center=no`, starting from the centering's shifts.
6. The lattice gate: `model_cavgs_rejection quality_model=nano` at the grade for raw class averages, 0.52.
7. The model gate, the second rejection, coupled with the first. The class averages the lattice gate kept are
   aligned to the start model (`refine3D_nano oritype=cls3D maxits=1`, a global search). A class is rejected when
   its correlation lies below `median - K_MODEL x 1.4826 x MAD` of the survivors (MAD, the median absolute
   deviation), and a class is kept only when both gates keep it. Then the prune, the original frame index recorded
   as `frame`, and the 2D shifts mapped to the 3D field.
8. The kept raw frames extracted.
9. The kept frames, only, flex-denoised.
10. The denoised frames swapped into the project.
11. `autorefine3D_nano` without the atomicity prior, `icm=no`, `lp=1.5`, with the centred mask.

`trajectory_solve3D.txt` records each stage's wall time, frames in and out, and numbers. A rerun with `mkdir=no` in
the execution directory resumes after the last finished stage. The model gate's override key is `model_gate_nsd`,
which the parameter layer refuses with any other program. Its products are the class table `model_gate_classes.tsv`
(class, population, lattice peak, both states, correlation, best orientation and shift) and class-average/
reprojection pairs.

## How the gates were set

**On denoised class averages (phase 0).** The maintainer's labelled classes of 71-17/traj04, 71-6/traj01 and
71-6/traj04 were rebuilt on the flex-denoised stacks: populations equal class by class, and the stacks bit-identical
to the previous run's. Each class was aligned to the simulated start model.

- The correlation ranks good from bad classes as well as the lattice peak does: area under the curve (AUC) 0.990,
  0.959 and 0.982, pooled 0.982, against 0.990, 0.975 and 0.983 for the lattice peak.
- The rule of ruling 2 (the smallest K in 1.5, 2, 2.5, 3, 4 at which, behind the lattice gate, at least 98 percent
  of the selected classes are kept on every set) gave `K_MODEL = 3`: recall 0.984, 0.988 and 0.985, while 2.5
  falls to 0.979 on 71-6/traj01.
- At 3 the cut removed 0, 9 and 0 of the lattice gate's survivors; the 9 were all rejected by the maintainer.
- Tables: `run/scratch/keep/gate/` (`classes.tsv` and `gates.tsv` per set, `gates_pooled.tsv`, `auc.tsv`); montages
  in `run/scratch/keep/gate/<set>/montages/`.

**On raw class averages (phase 1, ruling 6).** The same classes were rebuilt on the raw frames
(`run/scratch/keep/gate_raw/`, `lattice_scan.tsv`).

- A defect was found first and fixed: the nano validity gates judged flatness after a 10 Å low-pass, and so
  hard-rejected raw averages, including 61 classes the maintainer selected.
- After the fix the lattice peak still ranks the raw averages (pooled AUC 0.978).
- But 0.65 keeps only 0.984, 0.962 and 0.939 of the selected classes. The same rule as for 0.65 gives 0.52 (recall
  0.995, 0.995, 0.985; removal of the rejected 0.22, 0.20, 0.37), now `LAT_PEAK_MIN_RAW`, which the program passes
  as `thres`.
- Behind it every K keeps the same recall. The rule fell to the strictest, 1.5, which removed two of thirteen
  classes of a clean simulated trajectory. Ruling 9 broke the tie toward recall: `K_MODEL = 3`.

**What that means.** The model gate is nearly inert: at 98 percent recall, a cut relative to the survivors' median
and spread has no room. The correlation level to the simulated particle also differs between trajectories
(survivors' medians 0.26 to 0.58), so an absolute floor would not transfer. A gate relative to the upper mode of the
correlations is the form to try next.

## What gate 2 removed on the ten trajectories (phase 2)

Columns, in order: classes and frames the lattice gate removed; classes and frames the model gate removed; frames
kept, with the previous run's value; then the 3D numbers with the previous run's beside them (Fourier shell
correlation (FSC) 0.143 of the final halves at the program's mask, which differs from the previous run's by up to
5 Å).

| trajectory | lattice gate | model gate | frames kept [prev] | FSC 0.143, Å [prev] | strong atoms [prev] |
| --- | --- | --- | --- | --- | --- |
| 71-17/traj01 | 0 / 0 | 0 / 0 | 12,971 [12,846] | 1.12 [1.36] | 374 [364] |
| 71-17/traj02 | 10 / 129 | 0 / 0 | 12,842 [11,646] | 1.12 [1.12] | 324 [401] |
| 71-17/traj03 | 1 / 2 | 0 / 0 | 12,969 [12,708] | 1.12 [1.12] | 605 [549] |
| 71-17/traj04 | 31 / 445 | 6 / 133 | 12,393 [11,235] | 1.36 [1.36] | 222 [237] |
| 71-17/traj05 | 25 / 302 | 2 / 38 | 12,631 [9,536] | 1.40 [1.36] | 294 [298] |
| 71-17/traj06 | 8 / 87 | 0 / 0 | 12,884 [12,618] | 1.30 [1.30] | 655 [553] |
| 71-6/traj01 | 70 / 1,160 | 0 / 0 | 12,502 [8,186] | 1.13 [1.13] | 224 [209] |
| 71-6/traj03 | 116 / 2,148 | 2 / 35 | 11,479 [4,345] | 1.34 [1.65] | 98 [172] |
| 71-6/traj04 | 140 / 2,811 | 4 / 81 | 10,770 [5,890] | 1.23 [1.26] | 178 [195] |
| 71-6/traj02 | 8 / 101 | 0 / 0 | 13,561 [11,718] | 1.33 [1.10] | 212 [248] |

- **Gate 2 removed little.** It removed 14 classes and 287 frames in all, on five trajectories. Its class-average/
  reprojection pairs are in `run/scratch/keep/pipeline/<dataset>/<trajectory>/<trajectory>/1_trajectory_solve3D/7_refine3D_nano/model_gate_rejected_pairs.mrc`
  for 71-17/traj04, 71-17/traj05, 71-6/traj03 and 71-6/traj04, beside the class table and twenty kept pairs.
- **The lattice gate's montages** are `6_model_cavgs_rejection/quality_selected_cavgs.mrc` and
  `quality_rejected_cavgs.mrc` in the same directories. For 71-17/traj01 the directory is `2_trajectory_solve3D`.
- **Far fewer frames rejected.** Together the gates removed 5.6 percent of the frames, against the previous run's
  24 percent. On 71-6/traj03 and traj04 that means 2.6 and 1.8 times as many frames in 3D, with maps no worse.
- **The 3D otherwise mostly unchanged.** The 3D numbers equal the previous run's within 0.04 Å on seven
  trajectories. 71-6/traj02 is worse (1.33 against 1.10 Å, with a mask 1.8 Å larger), and 71-17/traj01 better
  (with a mask 5 Å smaller).
- **The shifts are tamed.** Centering first and refining with `center=no` keeps the largest shift at 26 to 37 px,
  against up to 83 px from the previous run's repeated recentering.
- **In situ on the labelled trajectories** (the workflow's own classes labelled by the majority of their frames'
  hand labels), the gates keep 98.5 to 100 percent of the mostly-selected classes and remove 6 to 25 percent of the
  mostly-rejected ones (`run/scratch/keep/pipeline/insitu.tsv`).
- The full table, with the previous run's values beside, is `run/scratch/keep/pipeline/table.tsv`.

## The remaining tail per frame

The plan's measure is the share of frames whose correlation in the 3D (`ptcls_vs_reprojs_corrs.csv`) lies below
`median - 3 x 1.4826 x MAD`, and the longest run of consecutive such frames. It is empty on nine trajectories and 4.9
percent on 71-17/traj03 (the longest run 23 frames). The per-frame correlations are broad (median 0.26 to 0.58, with a
spread comparable to the median), so this cut lies below nearly every frame. A frame-level gate after the 3D, which
ruling 4 names as the remedy for a tail inside accepted classes, needs a different measure than this one.

## Every fix

Each fix was made in its owning module, in the phase that found it.

- **`make_cavgs` writes `cls3D`** (phase 0, `simple_make_cavgs_strategy.f90`). It made the class-average field
  congruent in memory but wrote only `cls2D` and `out`, so `refine3D` on `cls3D` crashed on its project.
- **A failed run counted as passed** (phase 1, `simple_parameters_phases.f90`; and under ruling 4,
  `simple_parameters_registry.f90`). A missing input volume or stack, or a refused file format, stopped with a plain
  `stop`, exit status 0, so CTest counted the failure as a pass. These now stop with `THROW_HARD`.
- **`refine3D` on an empty `cls3D`** (phase 1, ruling 4, `simple_refine3D_strategy.f90`) now stops with a message
  instead of a segmentation fault.
- **The nano validity gate's flatness** (phase 1, `simple_cavg_quality_feats.f90`) is measured on the unfiltered
  class average, with a tester case.
- **`autorefine3D_nano` only defaults `mkdir=yes`** instead of forcing it (phase 1), so a workflow can run it in its
  own stage directory. Inside `simple_test_exec` no program makes an execution directory.
- **`model_cavgs_rejection` writes the raw lattice peak** on every scored `cls2D` row, for the class table
  (phase 1).

## What changed, phase by phase

The working tree held the previous run's two fixes and its plan note when this run started (the `make_cavgs` even/odd
halves, `autorefine3D_nano`'s transfer table, and its test assertion). They are in the complete diff but are not this
run's work.

- **Phase 0** (`run/review/phase_0.diff`): the `make_cavgs` fix and the phase 0 Progress entry.
- **Phase 1** (`run/review/phase_1.diff`):
  - the program (commander in `single_commanders_nano3D.f90`, UI entry in `single_ui_trajectory.f90`, router
    case);
  - the key `model_gate_nsd` and the program's parameter checks;
  - `model_gate_states`, `K_MODEL`, `LAT_PEAK_MIN_RAW` and the tester cases in `simple_cavg_quality_lattice*.f90`;
  - the fixes above;
  - the `single_workflow` test: `fcc` runs the program on the simulated trajectory, and `wurtzite` keeps the hand
    chain and now carries the atomicity-prior refinement `fcc` had;
  - the policy note `model_cavgs_rejection.md` and the test descriptions.
- **Phase 2** (`run/review/phase_2.diff`): the plan's Progress only; no source changed.
- **Phase 3**: this report and the plan's Progress entry.

## Decisions taken under delegation

Each is recorded in the plan's Progress, under the phase named.

- **Phase 0:**
  - the start model from the diameter of the labelled classes, with the program's chunk estimate as a check
    (both gave K = 3);
  - the montages as normalised pure-Python stacks.
- **Phase 1:**
  - stage directories named after the programs they run, records in `trajectory_solve3D.txt`, resume with
    `mkdir=no`;
  - diameters taken at the 0.1 Å `estimate_diam` reports;
  - `element` (default Pt) for the start model;
  - `maxits`, `maxits_between` and `nspace` forwarded to the 3D (the one-thread test needs them);
  - the model gate's alignment with `volrec=no`;
  - the atomicity-prior coverage moved to the `wurtzite` suite;
  - the original frame recorded as `frame` before the denoised stack is swapped in;
  - the lattice threshold on a 0.01 grid.
- **Phase 2:**
  - 71-17/traj01 not rerun, because its model gate removed nothing at 1.5 and so nothing at 3;
  - the keep rule (the last 2D iteration's class averages, the stacks dropped and remade on demand);
  - the per-frame tail reported with the median and 5th percentile beside it.

## Deviations from the plan, and why

- **The program's order.** It is rulings 6 and 7's raw-first order, not section 2's: denoising made junk look like
  particles.
- **There is no second diameter estimate or second centering** (ruling 7).
- **The lattice gate runs at 0.52 inside the program, not 0.65.** The rule that set 0.65 gives 0.52 on raw averages
  (ruling 6).
- **`K_MODEL` is 3 by ruling 9**, after the raw tables tied at every K.
- **The test has no planted junk** (ruling 8). The plan's noise windows, ruling 4's blobs and their assertions were
  tried and withdrawn. Planted junk is pooled into particles by the denoiser and mixed into particle classes by the
  2D, beyond what class-average gates can remove. The `fcc` suite asserts the map's floors, the pose checks and at
  most 5 percent of frames lost to the gates.
- **Phase 1 stopped three times**, and each stop was answered by a ruling (4, 5 to 8, 9). The records are in the
  plan's Progress, phase 1.

## Evidence

- **Tests** after the last source change:
  - the fast gate passes 15 of 15 (`run/scratch/keep/fast_gate_p1_7.log`);
  - `single_workflow_fcc` passes (`run/scratch/keep/ctest_fcc_p1_6.log`): no frame lost, map correlation 0.985
    with the truth, FSC 0.143 at 3.58 Å, pose checks within bounds;
  - `single_workflow_wurtzite` passes (`run/scratch/keep/ctest_wurtzite_p1_3.log`);
  - `unit_single`, `unit_ui`, the nano tester, and the registry and description checks pass.
- **Against the previous run:**
  - denoised stacks bit-identical (phase 0);
  - in the denoise-first form, the program's masks and lattice-gate states identical to the hand chain on
    71-17/traj01 (phase 1);
  - all ten trajectories tabled beside the previous run's values (phase 2).

## What is left open for the maintainer

- **Validate the ten maps** (`final_results/` of each `11_autorefine3D_nano`). 71-6/traj02's lower resolution
  (1.33 against 1.10 Å) is the one clear loss against the previous run.
- **The model gate is nearly inert.** The next form to try is a cut relative to the upper mode of the correlations.
- **A frame-level gate after the 3D,** for a tail inside accepted classes, needs a measure other than this run's
  median-relative cut, which finds almost nothing.
- **The flex denoiser turns structureless frames into pooled particles.** This is recorded, not acted on. It is why
  the program now gates before it denoises.
- **The fresh-start rule of `refine2D_nano`** leaves time-chunk classes only partly contiguous: 0.51 to 0.76 of
  consecutive frames share a class. Ruling 4 accepted it as a finding.
- **The identity outputs of `single_workflow`** in the test registry still name the hand chain's layout. The `fcc`
  suite's products now sit two levels deep (`4_trajectory_solve3D/11_autorefine3D_nano/`). Not changed in this run.
