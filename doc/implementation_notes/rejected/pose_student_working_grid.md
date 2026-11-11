# The pose student on the low-pass working grid

> **Negative result (2026-11-11).** The autonomous run `student_crop` stopped in Phase 0: the verdict rule's
> self-check (each of three seeds against the other two's range) cannot pass on any code. The change was then made
> interactively and measured (base six seeds, crop five to six seeds, `/data2/agent_runs/student_crop/scratch/keep/summary.tsv`):
> training ran 3.3 times faster on bgal, but the student lost on every fixed set (within 6° by 0.6 to 3 points; bgal's
> seeded fraction 0.417 to 0.284). It was not kept (patch: `/data2/agent_runs/student_crop/review/crop_interactive_2026-11-11.patch`).
> A cap on the network's angles alone at the same counts lost the same way, so the loss came from the angles. The
> kept change caps the network's half ring at twice π k_max (`doc/algorithms/pose_student.md`, "The network's half
> ring"). One premise below was wrong: the polar calculator's phase ramp uses the native box, so shifts are in
> native pixels at any crop.

Status, 2026-11-10: planned, for the autonomous run `student_crop` on the Dell. Base: the maintainer's commit of
the `student_speed` diff (with the product rule of 2026-11-10) and of this plan, on top of `13ad4833e`; the run
records the hash in Phase 0. Supersedes: none. Depends on: the pose student after `student_speed`
(`completed/pose_student_training_speed_report.md`, validated: yes, 3 seeds on 5 fixed sets) and solve3D's
low-pass crop rule (`lpstages_setlims` in `src/utils/filter/simple_estimate_ssnr.f90`, in production).

## 1. Why

The pose student always works at the native box. `learn_pose_student`, `predict_poses` and
`validate_pose_student` never set `box_crop`, so it defaults to `box`. The polar size `pftsz` follows the mask
radius at that box (`magic_pftsz(msk_crop)`, `simple_parameters_phases.f90`), whatever the low-pass limit. The
band `kfromto` already ends at `lp`, but the angles do not. The models trained on 2026-11-10 (seed 1, 24
threads) show this:

| set | box, sampling | `lp` | `pftsz` | `kfromto` | angles the outer ring needs (π k_max) |
| --- | --- | --- | --- | --- | --- |
| bgal | 256, 1.275 Å | 15 Å | 220 | 2-22 | 69 |
| medium assembly (T200 ... T) | 128, 3.0 Å | 15 Å | 160 | 2-26 | 82 |

The rest of SIMPLE works on a grid derived from `lp`. solve3D's stage schedule sets
`smpd_target = max(2 Å, lp/3)` and crops to the nearest magic box, at least 88 and never above the native box
(`lpstages_setlims` through `autoscale` in `src/utils/simple_magic_boxes.f90`). Since
`box * smpd = box_crop * smpd_crop`, a resolution keeps its Fourier index: `kfromto` is unchanged, and only the
angles and the image-side work shrink. At 15 Å the rule gives these grids:

| set | `box_crop` | `smpd_crop` | `pftsz` | factor on `pftsz` |
| --- | --- | --- | --- | --- |
| bgal | 88 | 3.71 Å | 72 (≥ 69) | 3.06 |
| medium assembly | 88 | 4.36 Å | 108 (≥ 82) | 1.48 |
| the CTest fixture (64 px at 2 Å, `lp` 8 Å) | 64 | 2.0 Å | unchanged | 1 (no crop) |

Two things scale with `pftsz`: the network's convolution positions and the readout's rotations (`nrots = 2 pftsz`).
So do the teacher's probability table, the polarisation and the prediction. The training tensors
(`pftsz x nk x 2 x nptcls`) shrink by the same factor. These are hypotheses that Phase 2 measures; the table is
the arithmetic of the rule, not a measurement.

## 2. Contract

- **Off state.** Where the rule gives `box_crop == box`, the programs behave byte for byte as on the base at a
  fixed seed: the model header, the weights file, the rows, the prediction table and the bench records apart
  from their timings. The CTest fixture is such a case (64 px at 2 Å with `lp` 8 Å: `smpd_target` 2.67 Å gives
  a box of 48, which the floors raise to 88, above the native 64, so `autoscale` keeps the native box).
- **One rule.** The student's crop is solve3D's rule, obtained by calling it, never by copying its constants.
- **The model carries its grid.** `box_crop` is already in the model header (`write_header`, `simple_pose_student.f90`
  `read_model`). `predict_poses` takes the crop from the model, so a model predicts on the grid it was trained on,
  and the existing geometry check (`student_geometry_mismatch`) keeps refusing a project whose grid differs.
- **Shift units.** The features, the model's `shift_sigma`, the augmentation and the prediction table work in
  native pixels or Å, as now. `polar_feats` builds its phase ramp as `2π/box` times the Fourier index, which is
  correct in native pixels on any crop because the index is shared. The polar calculator (`polarft_calc`) builds
  its ramp from its own image size (`argtransf`, `simple_polarft_core.f90`), so its shifts are in cropped pixels.
  Every shift that crosses between the two (the shift registration's start and result, any shift passed to
  `shift_ptcl` or the gradient search, the stored shifts the batch preparation applies) is converted exactly
  once.
- **What stays.** The architecture, the schedule, the readout, the HHCNN products, the teacher's search algorithm,
  the prediction table's columns and refine3D/solve3D are unchanged. The teacher's rows do change on a cropped
  grid, because the probabilistic table is filled on that grid; this is what refine3D at the same `lp`
  computes, and it is accepted (decision 7).

## 3. Decisions and rulings

1. (2026-11-10, Hans) The student works on the library's low-pass working grid. `learn_pose_student` derives
   `box_crop` from the project's `box` and `smpd` and the student's `lp` with solve3D's rule:
   `lpstages_setlims(box, 1, smpd, lp, lp, lpinfo, verbose=.false.)`, keeping `lpinfo(1)%box_crop`. That gives
   `smpd_target = max(2 Å, lp/3)`, a magic box, at least 88, never above `box`. It is called, not copied. If the
   call cannot be made from the student's layer (an import cycle, a side effect), stop and say why.
2. (2026-11-10, Hans) No new command-line key and no internal switch. `box_crop` stays off the UIs of
   `learn_pose_student`, `predict_poses` and `validate_pose_student`. The derivation sets it on the program's
   own command line before its single `params%new`, the way solve3D sets it on its stages' command lines.
3. `predict_poses` sets `box_crop` from the model header (`infile`) before `params%new`, never from its own
   derivation. A model of the base (header `box_crop == box`) predicts as before. The model format is unchanged,
   so the version stays 4.
4. `pftsz` follows `msk_crop` through `params` as everywhere else (`magic_pftsz`). The student adds no angular
   rule of its own.
5. The off state of section 2 is a gate: byte identity at a fixed seed where the rule gives `box_crop == box`.
6. (2026-11-10, Hans, as for `student_speed`) The student is judged on the fixed teachers' reads against the
   base's spread over three seeds, never on one pair of runs. The acceptance of equal reads is Hans's; the run
   sets no floor. The verdict rule of section 6 decides only whether Phase 3 runs.
7. The teacher's rows on the cropped grid are accepted as the teacher, because they are refine3D's at that `lp`.
   The truth sets judge the student against truth poses, so their reads stay comparable with the base's.
8. Shifts: native units outside the calculator, cropped pixels inside it, converted once at each crossing
   (section 2). A test plants a known shift and fails if a crossing's conversion is dropped or doubled
   (mutation-tested).
9. A rule other than solve3D's (a finer target such as `lp/4`) is measured only in Phase 3, only on its stated
   condition, and only as a variant outside the tree. Adopting one is Hans's ruling: the run stops with the
   table.

Defaults taken without a question: the run name `student_crop`; this plan's name; the five fixed sets, the seeds
1, 2 and 3 and the commands of `student_speed`; the model version unchanged; the Release build with
`-O3 -march=native` for all timings.

## 4. Not in this plan

- The network, the schedule, the design choice, the readout, the products of `simple_hhcnn.f90` (an edit there
  only if a cropped shape breaks a tile or block boundary; see the file table).
- The teacher's search (`eul_prob_tab`, the batch pipeline), refine3D, solve3D, `lpstages`' constants and
  `autoscale`.
- The prediction table's columns, the route, the student's consumers in refine3D (`prob_neigh_mode=geom`, the
  cones), `validate_pose_student`'s phases and reconstructions apart from the grid they run on.
- A 2D rule (`smpd_targets2D`) or any new grid rule.
- The Mac and the build without a team-safe BLAS are not timed; the HHCNN suite runs on the latter once
  (Phase 4), because the cropped shapes are new to its loops.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | base builds, fast gate, the base profile (3 seeds x 5 sets), the CTest entry, `validate_pose_student` twice |
| 1 | result-changing where the rule crops; identical where it does not | the unit suites, CTest `pose_student`, the off-state byte identity, the shift test with its mutation |
| 2 | measurement | the tree profile (3 seeds x 5 sets) against the base's spread; the verdict of section 6 |
| 3 | variants only (conditional) | stop with the table for Hans's ruling |
| 4 | end to end and close | `validate_pose_student` twice on tree and base, CTest, fast gate, `unit_numerics` on both BLAS builds, docs, report |

## 6. Acceptance

**The reads** are those of `student_speed` (its plan, section 1), from the same harness. On the truth sets
(`T200`, `T1000`, `T5000`, `T`, predicting the centred test half `S`): the folded mode within 6° of the truth
and within 12°, the median distance, the truth in the cone, the median cone size, the seeded fraction and its
in-cone share. On bgal: the training rows' agreement with the teacher (`dist_mode`, `in_cone`) and the seeded
fraction of the unseen rows. Also recorded: the held-out recovery and threshold from the model header, the
widths trained and kept, the epochs.

**The spread** (measured, Phase 0): for each set and read, the minimum and maximum of the base's three seeds.

**The verdict rule** (guessed, frozen here before Phase 0): a read is *inside* when the tree's three-seed mean lies
in the base's [min, max]. It is *at the edge* when it lies outside that range by no more than the range's own
width. Otherwise it is *outside*. Phase 3 runs when any read on any set is outside in the direction of a worse
student. Phase 0 checks that the rule is satisfiable by applying it to the base against itself (each seed
against the other two); if the base fails its own rule, stop.

**The speed** is reported, not gated. Recorded per set and seed: the wall time of `learn_pose_student`, its rows
fill, the designs' training split, peak memory; the wall time of `predict_poses` and its features-and-forward
time; the derived `box_crop`, `smpd_crop` and `pftsz`. The expected factors (guessed, from section 1's
arithmetic): training at least 2.5 times faster on bgal and at least 1.3 times on the medium sets. A set that
gains less is a finding to explain from the bench records, not a stop.

**Identity** (Phase 1): on the CTest fixture's geometry at a fixed seed, the base and the tree write identical
model headers, weights and prediction tables. The exclusions are the timing fields of the bench records and of
the header's time line, declared here.

**Reported, not gated:** `validate_pose_student`'s maps from the student's own orientations on bgal (tree and
base, two runs each) with their FSC against the control; the paths are in Progress.

## 7. Common rules and stop conditions

- No new parameter or key; no floor, tolerance or assertion relaxed; every new number is a named constant with
  its reason; a file outside the table gets a table row before the edit; fix as you go within the table.
- Variants are built from a copy of the tree under scratch, with their patches kept; only the kept change enters
  the diff.
- OpenMP: nothing new is expected. If a parallel region changes, the rules of `simple-modern-fortran` apply,
  with a 16-thread repeated identity check.
- Timings with the machine idle, `uptime`'s load recorded with each; a timing taken under load is repeated.

Stop when:
- decision 1's call cannot be made from the student's layer;
- the derivation needs a new key or a second code path;
- the off-state identity cannot be met;
- Phase 3's condition holds and its variants are measured (Hans rules on the rule);
- the base fails its own verdict rule in Phase 0;
- a unit suite of the base fails on the tree and the cause is outside this run's files;
- a file outside the table must change in a way the table's purpose does not cover.

Gates no ruling may waive: the fast gate; `unit_numerics` suites `hhcnn`, `pose_student`, `pose_student_strategy`
and `polar_feats`; CTest `pose_student`; the off-state byte identity; the shift test with its mutation; the
HHCNN suite's thread identity on both BLAS builds.

## 8. Phases

- **Phase 0, baseline.** Record the base hash, the machine and the load. Build the base export (Release) and the
  tree, run the fast gate, and build the variant without `SIMPLE_BLAS_TEAM_SAFE`
  (`student_speed`'s configure line). Copy the fixed sets' MD5 list and the harness from
  `/data2/agent_runs/student_speed/scratch/keep/` (sets read in place, never modified) and check the sums.
  Extend the harness's summary line with `box_crop`, `smpd_crop` and `pftsz` from the model header and the
  predict time split. Run the base profile: three seeds on every set. Run CTest `pose_student` with its wall
  time, and `validate_pose_student` twice on bgal (Phase 4's command). Compute what decision 1's rule gives for
  each set and the fixture, by hand from `lpstages_setlims`' formula, and write it into Progress; this is the
  table Phase 1's unit test must reproduce. Apply the verdict rule to the base against itself.
  Exit: `summary.tsv` holds 15 base lines with the new columns; the spread table and the rule's self-check are in
  Progress; the CTest time and the two validation runs are recorded.

- **Phase 1, the derivation.** Before editing, read how the polar calculator, the batch preparation
  (`build_batch_particles3D`), the reprojection model (`materialize`/`read_reprojection_model`) and the shift
  registration (`pftc_shsrch_grad`) handle `box_crop < box` in refine3D. Then trace every shift crossing in
  the student's strategy, and write both findings into Progress before the edit.
  - `learn_pose_student`: derive `box_crop` (decision 1) and set it on the command line before `params%new`;
    read the project's `box` and `smpd` the way existing commanders read them before `params%new` (search
    first).
  - `predict_poses`: take `box_crop` from the model header (decision 3).
  - `validate_pose_student`: nothing of its own; check that its learn and predict steps inherit the grid.
  - The learn bench record and the predict bench record print `box_crop`, `smpd_crop` and `pftsz` (they print
    `box / box_crop` already).
  - Tests:
    - the derivation on section 1's three geometries plus one large box (400 px at 1.0 Å, `lp` 15 Å),
      against Phase 0's table;
    - a model round trip with a crop;
    - `predict_poses`' refusal of a project whose derived grid differs from the model's, where that is
      testable;
    - the shift test of decision 8, with its mutation recorded;
    - the off-state identity at the fixture's geometry, base against tree.

  Exit: the fast gate; the four `unit_numerics` suites pass; CTest `pose_student` passes; the identity holds
  with the declared exclusions; the mutation fails the shift test and is removed; `T200` seed 1 trains and
  predicts on the cropped grid with the expected `box_crop`, `pftsz` and `kfromto` in its header.

- **Phase 2, the measurement.** Run the tree profile: three seeds on every set, with the harness. Fill the table
  of reads and times beside the base's, apply the verdict rule, and attribute the gain per set from the bench
  records (rows fill, training, prediction).
  Exit: 15 tree lines in `summary.tsv`; the comparison table and the verdict in Progress; Phase 3 marked "runs" or
  "skipped (no read outside)".

- **Phase 3, a finer target (only on its condition).** Build variants with `smpd_target = max(2 Å, lp/3.5)` and
  `max(2 Å, lp/4)` as patches under `scratch/variants/` (the student's derivation only). Run each on the sets
  whose reads were outside, three seeds each.
  Exit: the variants' table in Progress, then STOP for Hans's ruling (decision 9). When skipped, the exit is the
  Progress line that says why.

- **Phase 4, end to end and close.**
  - Run `validate_pose_student` twice with the tree and twice with the base: each in its own directory with a
    copy of `/data2/bgal_addon/3_selection/bgal_addon.simple`, using
    `prg=validate_pose_student projfile=bgal_addon.simple nsample=1500 mskdiam=180 pgrp=d2 pgrp_start=c1 nthr=24`.
    Record its phase timings, the readout reconstructions' FSC and the shift registration's records.
  - On the default build: CTest `pose_student` at its registered thread count, with its time. On both BLAS
    builds: the fast gate and `unit_numerics`.
  - Update `doc/algorithms/pose_student.md` with a paragraph on the working grid (the rule, where it is
    derived, that the model carries it).
  - Write the report as `doc/implementation_notes/completed/pose_student_working_grid_report.md`. Move this
    plan to `completed/` with `mv` and update its row in `doc/implementation_notes/README.md`.

  Exit: all of the above recorded in Progress; the report ends with a commit message.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/parallelization/simple_strategy_pose_student.f90` | the derivation, predict's grid from the model, shift crossings, bench records | 1, 3 (variants only), 4 |
| `src/main/strategies/parallelization/simple_strategy_pose_student_tester.f90` | the derivation and geometry tests | 1, 4 |
| `src/main/commanders/simple/simple_commanders_pose_student.f90` | validation's steps, if the grid must be passed on | 1, 4 (if needed) |
| `src/main/learn/simple_pose_student.f90` | model header and round trip, if needed | 1 (if needed) |
| `src/main/learn/simple_pose_student_tester.f90` | the round trip with a crop | 1 |
| `src/main/learn/simple_polar_feats.f90` | shift units, if a crossing is wrong | 1 (if needed) |
| `src/main/learn/simple_polar_feats_tester.f90` | the shift test | 1 |
| `src/main/learn/simple_hhcnn.f90` | only if a cropped shape breaks a tile or block boundary | 1, 2 (if needed) |
| `src/main/learn/simple_hhcnn_tester.f90` | a test of such a shape, if needed | 1, 2 (if needed) |
| `src/main/commanders/test/simple_commanders_test_highlevel_pose.f90` | the workflow test, if its checks need the grid | 1, 4 (if needed) |
| `src/main/commanders/test/simple_commanders_test_class.f90` | suite registration, if a tester is added | 1 (if needed) |
| `src/utils/filter/simple_estimate_ssnr.f90` | only to make `lpstages_setlims` callable without output, if it is not | 1 (if needed) |
| `doc/algorithms/pose_student.md` | the working-grid paragraph | 4 |
| `doc/implementation_notes/planned/pose_student_working_grid.md` | this plan, Progress | 0-4 |
| `doc/implementation_notes/completed/pose_student_working_grid.md` | this plan after the move | 4 |
| `doc/implementation_notes/completed/pose_student_working_grid_report.md` | the report | 4 |
| `doc/implementation_notes/README.md` | the row | 4 |

## Progress

(one entry per phase, written by the run)
