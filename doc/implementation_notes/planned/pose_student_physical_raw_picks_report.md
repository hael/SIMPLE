# Report: the pose student with the physical teacher and alignment over the cone

Run `pose_raw`, 2026-11-16 to 2026-11-18, on the Dell, unattended. Plan: `pose_student_physical_raw_picks.md` (beside
this report; its Progress section holds the detail). Bases: launched on `1ea0f8340`; moved to `fd9af0add` by ruling R4
(every measurement of Phases 0 to 4 was made there); moved to master's `68638be28` by ruling R6 (the final gates,
identity and validations). Outcome: **adopted under ruling R7's centring policy**, uncommitted for review.

The pose student is now trained on the physical posterior (no tempering), every row audited. Given the teacher's map,
`predict_poses` predicts each particle at its stored shift, then aligns it against the map's reprojections of every
direction of the predicted cone and its mirror, at every in-plane rotation over a 30 Å shift disk, and commits the
best alignment. The policy, in Hans's words: "centering comes from the solve2D run prior to solve3D on both particles
used for training and particles being predicted, shift refinement only enters after the prediction on the reference
obtained from the training set, thus mirroring the class student." Against today's route (prediction and a 5 by 5
shift registration at the 2D shifts) the unseen particles' map is as good or finer on every set, and on streptavidin
finer than the full search's own. The run's original idea, predicting on raw picks without 2D shifts, was measured and
not adopted. `validate_pose_student` judges against the physical posterior's best pose and reports the unseen map.

Terms. The pose student is a small network trained on a solve3D run (the teacher). For a new particle it gives a
posterior over 1,000 projection directions folded onto mirror pairs: its mode and its cone (the directions holding
99% of the mass, at most 64). The reference is the full search: each unseen particle's best pose under the physical
posterior over all directions, every rotation and the 30 Å disk, from the 2D shifts. The arms: C0, today's student
and registration at the 2D shifts (the yardstick); R0, the same on raw picks; S1 and SK, the alignment over the mode's
mirror pair or over the whole cone; A and P3, today's teacher (the search's tempering factor `prob_beta` times the
particle's temperature T) and the physical one (T alone); P3SK, the physical teacher on raw picks with SK; C3SK, the
same models at the 2D shifts with SK (the adopted route). The real sets are bgal, trpv1 and streptavidin, each with
three frozen solve3D teachers and three training seeds (nine pairs); the truth sets T200, T1000, T5000 and T are
simulated, with eight seeds. Comparisons are paired, with the paired standard error (SE) in brackets; "beyond 2 SE"
is the plan's rule. From ruling R5 on, the unseen map's Fourier shell correlation (FSC) decides: against the
teacher's map on the real sets, against the truth map on the truth sets.

## 1. What was learned

**Predict at the 2D shifts, then align over the cone** (Phase 4, ruling R7, `fd9af0add`). C0 → C3SK:

| C0 → C3SK | bgal | streptavidin | trpv1 |
| --- | --- | --- | --- |
| unseen map, FSC 0.143 (Å) | 5.12 → 4.32 (−0.80, 0.09) | 7.02 → 4.90 (−2.12, 0.38) | 7.74 → 7.77 (+0.02, 0.12) |
| unseen map, FSC 0.5 (Å) | 7.86 → 7.59 (−0.27, 0.07) | 9.51 → 8.58 (−0.93, 0.11) | 9.01 → 8.72 (−0.28, 0.09) |
| the full search's own map, FSC 0.143 (Å) | 4.31 | 5.06 | 7.78 |
| committed within 6° of the reference (%), diagnostic | 68.1 → 90.8 | 13.5 → 62.1 | 35.2 → 86.6 |

On the truth sets the map is finer at FSC 0.5 by 0.61 to 1.61 Å (T to T200), and the committed direction lies within 6°
of the truth for 97 to 99% of the test half (C0: 74 to 93%). The committed shift equals the reference's at the median
on bgal and trpv1: from the 2D shifts the alignment over the cone recovers the full search's answer.

**Raw picks lose what the 2D shifts give** (Phases 0 and 4). Today's student on raw picks (R0) holds the reference in
its cone 25 to 50 points less often than C0 on the real sets. With the physical teacher and SK on raw picks (P3SK) the
map was finer than C0's on bgal (−0.49 Å), equal on trpv1, but coarser on streptavidin by 1.91 Å at FSC 0.143 (SE 0.49),
so the raw-pick route was not adopted. The same models fed the 2D shifts (P3SK → C3SK) hold the reference in the cone
25 to 53 points more often, and streptavidin's map is 4.0 Å finer: the input shifts alone make the difference.

**Align over the cone (SK), not the mode's pair (S1)** (Phase 2, raw picks, today's teacher): SK commits the reference's
direction 6 to 50 points more often on every set and gives the finer map everywhere (0.14 to 2.79 Å at FSC 0.143).

**The physical teacher holds or improves the map while agreeing slightly less** (Phase 3, raw picks, SK):
streptavidin's map finer by 0.54 Å (SE 0.24), bgal and trpv1 within noise; the agreement reads dip by 0.5 to 1.6
points on four sets. This is the class student's pattern: a sharper teacher regularizes, and its maps do not suffer.
The audit on real rows found no row off the physical temperature (the largest mass difference 2e-4 of the bound).

**Costs** (reported, not gated; per particle, 24 threads, wall): the alignment over the cone 3.2 to 5.0 ms on the real
sets from the 2D shifts (3.4 to 6.2 from raw picks), beside the full search's table fill of 2.3 to 4.0 ms and C0's 0.08
to 0.39 ms. Its cost is set by the cone: 44% of bgal's cones, 90% of trpv1's and all of streptavidin's reach the
64-direction cap.

**The truth sets' FSC reads saturate**: at 3 Å per pixel FSC 0.143 crosses at Nyquist (6.00 Å) in every run of the
aligned arms, and C3SK's 0.5 crossing sits at 6.86 Å in every run. A mean FSC over the shells finer than 20 Å, added
as a diagnostic, still separates the arms.

## 2. What changed, phase by phase

The driver's diffs: `review/phase_0.diff` to `review/phase_4.diff` under `/home/elmlundho/agent_runs/pose_raw/`.

- **Phase 0** (2 files: this plan, `doc/implementation_notes/README.md`): the baseline on `fd9af0add` after ruling R4
  (the first, on `1ea0f8340`, kept as superseded): base builds, the off-state identity base against base, the
  uncentred test half, nine solve3D teachers built and frozen, C0 and R0 on every set, decision 9.
- **Phase 1** (13 files, two new): `simple_pose_student_align.f90` (new), the alignment after prediction over a list of
  directions with their mirrors through the existing polar shift map, with its tester; `predict_poses` aligning given
  `vol1` with `center=no`; the teacher's audit (`physical_row_audit`) and its test; the reference pass
  (`run_pose_posterior_mode`); `validate_pose_student` rewritten (frozen teacher through the existing key
  `dir_refine`, frozen references, decision 7's reads, the maps, the costs); the UI help; the unit suite's
  registration; `doc/algorithms/pose_student.md`; the distributed master's sigma2 bootstrap made identical to a
  shared-memory run.
- **Phase 2** (6 files): S1 removed (the strategy, the commander, the prediction table, the times file, the report).
- **Phase 3** (6 files): the tempered teacher removed (`learn_pose_student` always sets `prob_beta` 1, the audit always
  enforces it); decision 8 revised by ruling R5.
- **Phase 4** (the 13 files re-applied on `68638be28` by ruling R6, then 9 files): the centring policy of ruling R7.
  `validate_pose_student` predicts its unseen half at the 2D shifts; the dense search on the raw picks (D0) is gone
  with its frozen file and constant; the readout's shift read uses the stored shift; comments, the UI help and
  `doc/algorithms/pose_student.md` state the policy; this plan (decisions 3 and 8), the README row, this report.

## 3. Decisions taken under delegation

- How the reference is computed (decision 2): the dense rows at prob_beta 1 over the 30 Å disk from the 2D shift,
  because a raw-pick start put the mode within 6° of the truth for 72.5% instead of 99.1% (Phase 0).
- The point groups and start groups read from the projects' own solve3D jobs (ruling R2; Phase 0).
- `predict_poses` sets `center=no` (Phase 0: centring re-centred the map and moved every stored shift; Phase 1).
- The frozen teachers enter through the existing key `dir_refine`, its help in `simple_ui_solve3D.f90` (a file-table
  row added under ruling R3; Phase 1).
- Decision 8 (a): SK (Phase 2). Decision 8 (b): P3, on the maps (ruling R5; Phase 3). Decision 8 (c): the raw-pick
  route not adopted; under ruling R7, C3SK adopted (Phase 4).
- C3SK's real-set runs used a nine-line measurement variant of `fd9af0add`'s final tree (`variants/C3`): the unseen
  particles at their 2D shifts, the change the policy itself makes, and `infile` to predict with Phase 3's trained
  model instead of training (ruling R7: no retraining). The truth sets ran on `fd9af0add`'s build itself.
- The readout arm's shift read compares the stored shift, now the 2D shift (it compared zero, right for raw picks only).
- A validation that died of a segmentation fault in the base's `oris%get_pinds` parallel loop (a heap fault, as Phase
  0's malloc abort in solve3D) was rerun (Phase 2); master's `68638be28` links the Fortran runtime statically against it.
- The comparison tool leaves out undefined reads (no seeded particle), and a mean-FSC diagnostic was added where the
  truth sets' FSC reads saturate; neither decides (Phase 3).

## 4. Deviations from the plan

- **Three base moves and three rulings that changed the route.** Ruling R4 moved the base after Phase 1 (Phases 0 and 1
  redone, the identity on 6vxx only); ruling R5 made the unseen map the judge (a removal of P3 drafted on the agreement
  reads was undone unbuilt); ruling R6 moved the copy to master's HEAD before Phase 4 (it reached the rulings a minute
  after Phase 3's last session read them, so Phase 3 ended DONE and Phase 4's first session stopped instead); ruling R7
  withdrew decision 3 (no pre-centring) and replaced it with the centring policy after Phase 4 stopped at decision 8 (c).
- **The plan's premise did not hold**: dropping the 2D shifts was the point of the run; the alignment over the cone did
  not recover them on streptavidin. The physical teacher's gain in agreement seen on the class student did not carry
  over (its maps did).
- **Files outside the original table**: `src/main/ui/simple/simple_ui_solve3D.f90` (Phases 1 and 4, ruling R3); the
  algorithm document's row extended to Phases 2 and 3 (ruling R3). Rows left unused: `simple_polar_student.f90`,
  `simple_eul_prob_tab.f90`, the workflow test `simple_commanders_test_highlevel_pose.f90`,
  `doc/policies/test_environment_policy.md`, `polar_student_future_directions.md`.
- **Machine sharing** during Phase 2 (another session's jobs, load 30 to 55 for about two hours) inflated those walls.

## 5. Evidence against the baseline

- Off-state identity, `simulated_workflow` 6vxx at one thread, `SIMPLE_SEED=7`, `test=identity`: on `fd9af0add` base
  against base empty and the tree against the base 186 files the same, 0 differing after Phases 1, 2 and 3; on
  `68638be28` the final tree 185 and 0 (master's workflow writes one file fewer).
- The class student byte-identical to its base (model header, weights' checksum) after every phase and on `68638be28`.
- Fast gate 15 of 15, `unit_numerics` (2,423 assertions) and `unit_pftc_align2D3D` (6,324 on `fd9af0add`, 6,327 on
  `68638be28`) on the Release and Debug trees after every phase and once more after the last edit; CTest
  `pose_student` and `class_student` passed.
- Nine planted violations (the mirror not searched, the first reference kept, a window for the disk, the refinement
  ignored, a reference searched twice, a thread-dependent score, the audit's targets tempered, the masses unnormalised,
  a bound blind to the temperature) each caught by its suite, after every phase.
- `validate_pose_student` on the three real sets at teacher 1, seed 1 with the final tree on `68638be28`: unseen maps
  4.41 (bgal), 4.57 (streptavidin), 7.78 Å (trpv1) at FSC 0.143, beside the full search's 4.35, 4.16 and 7.78 Å and
  C0's 5.26, 8.58 and 7.98 Å at that pair on `fd9af0add`.
- Tables and logs: `scratch/keep/p2/`, `p3/`, `p4/`, `reads_real.tsv`, `reads_truth.tsv`; `scratch/logs/`.

## 6. Left open for the maintainer

- **Reported, not gated**: the alignment over the cone costs about as much as the full search's table fill (3.2 to
  5.0 ms against 2.3 to 4.0 ms per particle), because 44 to 100% of the real-set cones reach the 64-direction cap; the
  student seeds few real-set particles (0 to 3%). A narrower cone (a smaller cap or a tighter mass) would cut the cost; untested.
- **Data-dependent**: the truth sets' FSC reads saturate at 3 Å per pixel; a finer read (a smaller pixel, or a
  band-averaged FSC such as the diagnostic used here) is needed to separate close arms. One validation per real set on
  `68638be28` is not paired with `fd9af0add`'s: the base's training numerics moved between the two.
- **Follow-ups**: the stream, where the unseen particles' 2D shifts come from the class student's own alignment (the
  policy assumes it; not run here); whether the alignment over the cone could replace the neighbourhood search
  (`prob_neigh geom`) that consumes the cones in refine3D (the contract kept that route unchanged; untested).
- **Skill text proposed** (not applied): in `simple-main-strategies-parallelization`, "The pose student's
  `predict_poses`, given `vol1`, predicts each particle at its stored (prior solve2D) shift and aligns it after the
  prediction over the cone's directions and their mirrors (`simple_pose_student_align`, the polar shift map unchanged),
  committing the best alignment; `learn_pose_student` reads the physical posterior (`prob_beta` 1) and audits every row
  against it."

Changed and new files (`git status --porcelain`; the fast gate passed once more on both trees after the last edit,
`scratch/logs/p4_final_gates.log`):

```
 M doc/algorithms/pose_student.md
 M doc/implementation_notes/README.md
 M doc/implementation_notes/planned/pose_student_physical_raw_picks.md
 M src/defs/simple_defs_hhcnn.f90
 M src/main/commanders/simple/simple_commanders_pose_student.f90
 M src/main/commanders/test/simple_commanders_test_class.f90
 M src/main/learn/simple_pose_student.f90
 M src/main/learn/simple_pose_student_tester.f90
 M src/main/strategies/parallelization/simple_strategy_pose_student.f90
 M src/main/ui/simple/simple_ui_solve3D.f90
 M src/main/ui/simple_test/simple_test_ui_class.f90
?? doc/implementation_notes/planned/pose_student_physical_raw_picks_report.md
?? src/main/strategies/parallelization/simple_pose_student_align.f90
?? src/main/strategies/parallelization/simple_pose_student_align_tester.f90
```

## 7. Commit message

```
pose student: physical teacher, alignment over the cone after prediction

The pose student follows the class student's centring policy: training
and predicted particles are centred by the prior solve2D, and shift
refinement enters only after the prediction, over the cone, against the
training set's map.

learn_pose_student trains on the physical posterior (prob_beta 1) and
audits every row against exp(-T (L - Lmin)) within the rounding bound of
the single-precision losses. Given the teacher's map, predict_poses takes
each particle at its stored shift (center=no), predicts the cone and
aligns the particle against the reprojections of every cone direction and
its mirror, every in-plane rotation over a 30 A disk, through the existing
polar shift map (simple_pose_student_align, new, with its tester); the
best alignment is committed with the dynamic centring's state. The 5x5
shift registration is gone. validate_pose_student takes a frozen teacher
(dir_refine), predicts its unseen half at the 2D shifts, judges against
the physical posterior's best pose and reports the unseen maps and the
costs. The distributed master's sigma2 bootstrap now matches a
shared-memory run.

Measured (doc/implementation_notes/planned/pose_student_physical_raw_picks.md):
the unseen maps are as good or finer than today's route on every real and
simulated set (streptavidin 7.02 -> 4.90 A at FSC 0.143); the alignment
over the cone beats the mode's pair; the physical teacher's maps are as
good or finer; predicting on raw picks was measured and not adopted.
Off-state identity and the class student unchanged.
```
