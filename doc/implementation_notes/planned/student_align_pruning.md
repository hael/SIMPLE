# The students' alignment: refine the screen's best, not every candidate

Status, 2026-11-18: planned, to run on the dell under arun (run `align_prune`). Base: `a26f57437` with the working
tree's per-thread feature buffers (pinned at launch; byte-identical to the base, checked by hand 2026-11-18).
Supersedes: none. Depends on: the class student on raw picks with the alignment over its candidate list (validated:
yes, `completed/class_student_raw_picks.md`); the pose student's alignment over the cone at the 2D shifts (validated:
yes, ruling R7 of `planned/pose_student_physical_raw_picks.md`).

## In one paragraph

After prediction both students align each particle against every candidate: the class student over its list of up to
16 classes, the pose student over its cone of up to 64 directions and their mirrors. Every candidate gets a coarse
shift map over the disk and a continuous refinement, because both call the polar shift-map search with `nbest` equal to
the list's length. The classic 2D search calls the same search with a small `nbest`: it screens a pool of
max(2·`nbest`, the screen's posterior-mass size) candidates, ranked by their score at the stored shift, and refines only
the screen's best `nbest`. The alignment is the students' main cost (class 0.5 to 1.2 ms, pose 3.2 to 5.0 ms per
particle, about the full search's own table fill). This run measures what a smaller `nbest` costs in class averages and
unseen maps, and adopts the smallest that loses nothing.

## 1. Why

- The costs above (`completed/class_student_raw_picks_report.md`; `planned/pose_student_physical_raw_picks_report.md`).
  Most pose cones reach the 64-direction cap (44 to 100% of particles on the real sets), so the pose student aligns
  against up to 128 references per particle.
- Pruning ranks candidates by their score at the stored shift. For the pose student under R7 the stored shift is the 2D
  solution's, so that score is informative. For the class student on raw picks it is weaker (a particle 7 to 9 Å off
  centre), which is why aligning over the whole list beat the most likely class alone; pruning may hurt there.
- Maps and class averages decide (ruling R5 of the pose run); agreement reads are diagnostics.

## 2. Contract

- Phases 0 to 2 change nothing in the tree's source: the arms are variant patches on `git archive` exports under
  `scratch/keep/variants/` (a developer key in the variant only). Phase 3 changes the tree only if a verdict adopts a
  smaller `nbest`, and only as section 8's file table says.
- No new command-line key in the tree: an adopted value is a named constant in `src/defs/simple_defs_hhcnn.f90`.
- The shift-map search (`simple_shiftmap_srch`) is read, not edited; the classic searches stay byte-identical.

## 3. Decisions

1. **Data and inputs.**
   - Class student: the three import projects of the class runs, at their settings
     (`/data/10335-streptav/tst/1_import_particles/tst.simple` `nsample=5000 ncls=30 mskdiam=85`;
     `/data2/bgal_addon/1_import_particles/bgal_addon.simple` `2500 50 180`;
     `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple` `2500 50 175`), eight seeds per set.
   - Pose student: the pose run's frozen teachers and physical-teacher models (`/data2/agent_runs/pose_raw`, read only;
     nine teacher-and-seed pairs per set), predicting the unseen half at its solve2D shifts as that run's arm C3SK did,
     on the variant sources of `/data2/agent_runs/center_ceiling/scratch/keep/variants/V2` (the C3 sources linked with
     `-static-libgfortran`; copy them).
   Checksums recorded in Phase 0; a changed or missing input is a stop.
2. **Arms.** `nbest` is the number of candidates the shift-map search refines (`smap%srch`'s third argument):
   - class student: all (today), 4, 2, 1, run inside one `validate_class_student` per set and seed (one teacher,
     model, reference and unseen set; the arms after today's, as the center_ceiling run's variant V1 did);
   - pose student: all (today), 16, 8, 4, inside one `validate_pose_student` per pair.
   Every arm records, per particle, the candidates screened and refined.
3. **Reads.**
   - Class student: the class averages' interpolated FRC 0.5 crossing, the median over the classes compared in every
     arm (at least `CLASS_MINPOP_FRC` unseen members in the reference and in every arm; the center_ceiling run's
     ruling R4), the mean FRC to 15 Å beside it; diagnostics: the committed class's agreement with the reference, the
     committed shift's error; the alignment's cost (ms per particle, wall, 24 threads).
   - Pose student: the unseen map's FSC 0.143 crossing (Å) against the teacher's map, FSC 0.5 beside it; diagnostics:
     the committed direction within 6° of the reference, the committed shift's error; the alignment's cost.
4. **Verdicts, written now.** One-sided 95% bounds of the paired mean difference against today's arm: eight seeds for
   the class student (t₇ = 1.895), nine pairs for the pose student (t₈ = 1.860). The loss L = read(arm) − read(all) in
   Å (positive: the pruned arm is worse). Margin δ = 0.3 Å for both.
   - Per student, an arm is *non-inferior* if UB(L) < δ on every set.
   - The adopted `nbest` is the smallest non-inferior arm whose alignment cost is at most 0.7 of today's on every set
     (a gain worth a constant); otherwise today's stays and the run records the table.
   - The two students are decided separately.

## 4. Not in this plan

- The shift-map search itself, the disk radius, the candidate lists' size (`KMAX_CLASS`, the cone's cap) or mass.
- The stream link (the pose student at the class student's aligned shifts): a phase of its own if the maintainer adds
  it by a ruling.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement setup | base Release build and fast gate; input checksums; the class and pose variants built (fast gate each); each variant's `nbest`=all arm byte-equal to the unpatched program (class: one set at seed 1, one thread; pose: one pair against its recorded C3SK) |
| 1 | measurement | the class student's 24 validations; its verdict |
| 2 | measurement | the pose student's 27 pairs; its verdict |
| 3 | adoption or record | if a verdict adopts: the constants in the tree, the testers, the workflow tests `class_student` and `pose_student`, fast gate, `unit_numerics`, `unit_pftc_align2D3D`, refine2D and refine3D byte-identical (`simulated_workflow` 6vxx, `test=identity`, one thread), the docs; the report either way |

## 6. Acceptance

- Floors: none carried; every compared number is measured in this run on one build per student. The margin and the
  cost ratio are written here (guessed) and never moved without a ruling.
- Reported, not gated: the screened and refined counts, the per-step costs.

## 7. Common rules and stop conditions

- The thread and disk policies of `doc/how2s/how2run_autonomous_agents.md`. The maintainer works on the machine: never
  his checkout or its `build/`; record `uptime`'s load with every timing.
- **Stops:** an input missing or changed; a Phase 0 equality check failing; an edit needed in the shift-map search or
  outside the file table; refine2D or refine3D not byte-identical in Phase 3.

## 8. Phases

- **Phase 0, the setup.** Pin the base, build Release, fast gate; checksums; the two variants with their equality
  checks. Exit: the variants and checks in Progress.
- **Phase 1, the class student.** Three sets × eight seeds, four arms each; paired tables; the verdict.
- **Phase 2, the pose student.** Three sets × nine pairs, four arms each; paired tables; the verdict.
- **Phase 3, adoption or record.** For an adopted `nbest`: `CLASS_STUDENT_ALIGN_NBEST` and/or
  `POSE_STUDENT_ALIGN_NBEST` in `simple_defs_hhcnn.f90`, passed to `smap%srch` by the students' aligners; their
  testers check the pruned path; section 5's gates; `doc/algorithms/class_student.md` and `pose_student.md`. The report
  `student_align_pruning_report.md` beside this plan, which stays in `planned/`, with `git status --porcelain`.

### File table

| File | What | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/student_align_pruning.md`, `student_align_pruning_report.md` (new) | this plan, Progress, report | 0 to 3 |
| `doc/implementation_notes/README.md` | the plan's row | 0, 3 |
| `src/defs/simple_defs_hhcnn.f90` | the adopted constants | 3 |
| `src/main/strategies/parallelization/simple_class_student_align.f90`, `simple_class_student_align_tester.f90` | the class student's `nbest` | 3 |
| `src/main/strategies/parallelization/simple_pose_student_align.f90`, `simple_pose_student_align_tester.f90` | the pose student's `nbest` | 3 |
| `doc/algorithms/class_student.md`, `doc/algorithms/pose_student.md` | the alignment's description | 3 |

## Progress

(one entry per phase, written by the run)
