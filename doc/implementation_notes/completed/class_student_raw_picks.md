# The class student on raw picks: no centring, alignment after prediction

Archived in `completed/` on 2026-11-17 by maintainer decision. The status and open items below are as
last written.

Status, 2026-11-16: adopted (ruling R4), uncommitted for review. The class student predicts on raw picks and aligns
each unseen particle over the shift disk against every class of its candidate list; `center_ptcls` left the route
and the validation; the amplitude channel was rejected (decision 7 (b)). The report is
`class_student_raw_picks_report.md`. Base:
`1637eeb1c` (the class student learns the physical posterior, P of `student_inplane_posterior.md`). Supersedes: decision 9 of `completed/class_student_2D.md` as revised on
2026-11-13 (unseen particles centred by `center_ptcls` before the student labels them). Depends on: the class
student (`doc/algorithms/class_student.md`), `center_ptcls` (kept as a program, dropped from the route).

## 1. Why

The class student exists to extend a 2D solution to new particles without a 2D analysis of them. Since
2026-11-13 its route first centres every unseen particle with `center_ptcls`, a 2D analysis of its own: on the
maintainer's runs of 2026-11-16 the centring took 6.9 of 7.1 ms per particle on streptavidin, 1.6 of 2.1 on bgal,
1.2 of 1.6 on trpv1, more than the full search it was meant to replace. The premise recorded then ("mirror the pose
student, which always worked on centred particles") misdescribed the pose student: it predicts first and registers
the shift afterwards against the reprojection at its predicted pose (its validation inputs carried a 2D solution's
shifts; its pipeline has no centring step). The maintainer's feedback on all three sets: the student's class
averages are "in desperate need of shift refinement"; today's registration is a local gradient refinement of the
shift at the predicted class and angle (`predict_classes_exec`, `predict_one`), not a search over the picking error.

The maintainer's design (2026-11-16): no centring; the student predicts on the raw picks; then a single-reference
alignment per particle against the training chunk's class average of its predicted class, every in-plane rotation
and shifts over the disk. That turns the shift problem from a solve2D run into one small alignment per particle.
Searching each class of the student's list instead of the predicted class alone uses the network as a
neighbourhood structure (Hans, 2026-11-16: keep this in mind for the pose student's post-prediction refinement,
`planned/polar_student_future_directions.md`). Separately, the amplitude channel, shift-invariant by
construction, may let the network rank classes on uncentred particles at all (Phase 3: it did not, decision 7 (b)).

What was measured on 2026-11-13 (`completed/class_student_2D.md`, the shift probes) and frames the risk: picking
errors of 8 to 10 Å median (16 to 22 Å at the 90th percentile); the network's mode degrades at 5 to 10 px
off-centre; wider shift augmentation costs accuracy; power spectra alone learn nothing. Whether the student's
candidate list still holds the right class on raw picks, so that a wide alignment over the list recovers it, was
never measured: that is Phase 0's first question.

## 2. Contract

- Off state: refine2D and refine3D byte-identical at one thread with a fixed seed (`simulated_workflow`, both
  suites, `test=identity`); the pose student byte-identical to the base (model header and weights).
- `center_ptcls` stays as a program, unchanged; only the class student's route and its validation stop calling it.
- The training set stays the teacher's (solve2D's) particles at the teacher's shifts; only unseen particles change.
- No new command-line key or parameter; no second code path kept for a rejected variant. The amplitude variants
  live as patches on an export of the base until their verdict.
- The networks keep their symmetries (rotation equivariance on the grid, the conjugation symmetry; the amplitude
  path cannot tell r from r + 180°, the complex path does) and training stays bit-identical at a fixed thread count.
- Costs are reported, measured, and enter verdict 7 (c) only as stated there; small data sets close nothing
  (Hans, 2026-11-15).

## 3. Decisions

1. (2026-11-16, Hans) No pre-centring of unseen particles in the class student's route. They enter at their
   picked positions: their shifts set to zero, as imported.
2. (2026-11-16, Hans) **Alignment after prediction**, against the teacher's class averages (the training chunk's),
   every in-plane rotation of the polar grid and shifts over the disk of the validation's reference pass
   (`CLASS_DENSE_SHIFT_SIGMA`, three sigma), at the student's band. Two variants, both measured (Hans: "measure
   both"):
   - **S1**: the predicted class (the mode) alone;
   - **SK**: every class of the student's candidate list (its `ncand` classes, at most `KMAX_CLASS`), the best
     alignment kept, its class committed.
   Every unseen particle is aligned (shortlisted or not), so the reads cover both routes; the route field still
   says whether the particle was shortlisted. The alignment reuses the search code the 2D table's fill uses for
   the disk (the shift map and the gradient shift search of `src/main/pftc`), unchanged; a new routine there is a
   stop (section 7).
3. **The amplitude channel for the class student**: the design then planned for the pose student (the modulus of
   the complex input after the CTF multiplication and the shell scaling, its own three blocks of the same shape with
   a ReLU, its own readout of even angular orders against a per-class amplitude template, scores added before the
   maximum, or before the log-sum-exp of P's readout). Variants: B (P plus the amplitude path), B' (B with
   each particle's amplitude channel from another particle, a fixed permutation: the capacity control), C (the
   amplitude path alone: a wiring check that must reproduce the old negative, held-out classes recovered at most
   0.10 on bgal; guessed).
   (Revised in Phase 3 under ruling R2, a deviation for the maintainer to confirm: C recovered 0.484 of the held-out
   classes on bgal at seed 1, against the guessed 0.10. The old negative was the pose student's within-6° read, which
   needs the in-plane angle that amplitudes cannot give; a 2D class is recognisable from its amplitude spectrum. The
   patch shows exactly the path's contract instead, so the wiring check became that contract, its thresholds set
   after C's first run: C's predictions are shift-invariant (its reads at the picked positions and at the
   reference's shifts within one point of each other; measured 0.2, 0.1 and 0.0) and its read-out angle cannot
   resolve the half turn (of the held-out rows of agreeing class, at least four in five within 10° of the teacher's
   angle or of its opposite, and at least one in four at the opposite; measured 0.896 and 0.563, B 0.025 at the
   opposite).) Confirmed by Hans (ruling R4, 2026-11-16): amplitudes alone recognise 2D classes; the check is the
   design's properties, the shift-invariant reads and the half-turn ambiguity.
4. **Data**, the maintainer's settings of 2026-11-16, raw picks of each import project, read only (copy, never write):
   - streptavidin: `/data/10335-streptav/tst/1_import_particles/tst.simple`, `nsample=5000 ncls=30 mskdiam=85`;
     the maintainer validated on 11,402 particles (5,000 + 6,402): if the import project's active count differs,
     take `/data/10335-streptav/tst/3_selection/tst.simple` with its shifts zeroed and record which (Phase 0: the
     import project holds 11,402 particles, all selected, all at zero shift: the maintainer's count, taken);
   - bgal: `/data2/bgal_addon/1_import_particles/bgal_addon.simple`, `nsample=2500 ncls=50 mskdiam=180`;
   - trpv1: `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple`, `nsample=2500 ncls=50 mskdiam=175`.
   `validate_class_student` makes its own teacher per seed (as R6 of the last run did), the unseen particles at
   zero shift.
5. **Arms.** C0: the base route (P, `center_ptcls`, today's local registration), the yardstick. R0: P on raw picks
   with today's local registration (what dropping centring alone does). RS1, RSK: P on raw picks with alignment S1,
   SK. B, B', C: the amplitude variants on raw picks with the alignment Phase 2 chooses.
6. **Reads**, per set and seed, on the unseen particles, against the reference (the physical posterior's best
   class, from the dense 2D table with its shift disk on the same raw picks, as now):
   - particles given a shortlist; shortlists holding the reference's class; median list length;
   - the candidate list holding the reference's class (all particles), the mode's agreement with it;
   - after the alignment: the committed class's agreement with the reference, the shift's distance to the
     reference's shift (median and 90th percentile, Å), the in-plane angle's distance;
   - class averages of the committed assignment (every unseen particle at its aligned class, angle and shift)
     against the teacher's, FRC 0.5 interpolated and the band's mean FRC to 15 Å (as
     `completed/class_student_vs_reference_matching.md`), beside the reference's own averages;
   - cost per particle (ms, 24 threads, the machine's load recorded): prediction, alignment, and their sum against
     the reference's full search.
7. **Verdicts, written now** (paired seeds 1 to 8, the paired mean of the difference and its standard error; a read
   is worse when the difference lies in the worse direction beyond twice its standard error):
   - (a) **S1 or SK**: the one with the higher committed-class agreement; within 2 SE, the cheaper (S1).
     (Phase 2: SK, its aligned class agreeing with the reference's 7 to 38 points more often on every set.)
   - (b) **The amplitude channel**: B is adopted unless it is worse than the chosen raw-pick arm on a read of
     decision 6, and only if B' does not show B's gains (a gain B' also shows beyond 2 SE is capacity, and B is
     rejected); C must pass its wiring check or the patch is wrong (a stop until fixed).
     (Phase 3: B rejected, worse than the raw-pick arm on reads of every set, its bgal gain shown by B' too.)
   - (c) **Dropping centring** from the route and the validation: the best raw-pick arm (chosen alignment, with or
     without B) against C0, on every set: committed-class agreement and class-average FRC not worse than C0's, and
     its cost per particle (prediction plus alignment) below the reference's full search. Otherwise stop with the
     table for the maintainer (the raw-pick route is then not adopted).
     (Phase 4: not met as written. Agreement and class averages better than C0's on every set; the cost below the
     full search on bgal and trpv1, above it on streptavidin, 0.537 against 0.442 ms. Ruling R4 (Hans, 2026-11-16)
     waives the cost condition for streptavidin, a 30-class set where the full search is already cheap (0.44 ms) and
     small data sets decide nothing, and adopts the raw-pick route as built: no centring, the alignment after
     prediction against the candidate list (SK), no amplitude path. A cap on the list relative to the class count is
     a possible follow-up, not done.)
   Margins: twice the paired standard error (measured); the C wiring bound 0.10 (guessed). (Phase 3: the C bound replaced, see decision 3.)
8. **Early stop**: after Phase 0, if R0's candidate list holds the reference's class for fewer than half of the
   particles on every set, no alignment over the list can recover the class: record it, skip RSK, and continue with
   RS1 and the amplitude variants (the amplitude channel is then the only route left to the network).
   (Phase 0: not met on any set, R0's list holds the reference's class for 89.5 to 96.4% on average and for at
   least 88% in every run; RSK runs.)

Defaults taken without a question: eight paired seeds at 24 threads; S1 and SK measured in the same validation run
(the alignment runs twice on the same prediction); the reads' FRC tools of the reference-matching record.

## 4. Not in this plan

- The pose student (the neighbourhood idea goes to
  `planned/polar_student_future_directions.md` as a note only).
- `center_ptcls` itself, solve2D, refine2D, the stream (`planned/class_student_stream.md`).
- Any file outside the table.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | C0 and R0 at eight seeds on three sets, frozen |
| 1 | result-changing for the class student's prediction and validation | testers of the alignment (a planted shift recovered over the disk; SK keeps the best class); plants; fast gate, `unit_numerics`, `unit_pftc_align2D3D`; `class_student` workflow test; off-state identity |
| 2 | experiment | decision 7 (a); decision 8 if it applies |
| 3 | patch outside the tree, then experiment | C's wiring check; the HHCNN tester's symmetries and thread identity on the amplitude path; decision 7 (b) |
| 4 | adoption or record | decision 7 (c); section 5's Phase 1 gates on the adopted tree; `validate_class_student` on the three sets |

## 6. Acceptance

- Spreads and the paired standard errors come from the runs (measured); the C wiring bound is guessed and says so.
- Off-state identity: `simulated_workflow` both suites on the base and on the final tree at one thread, seed fixed,
  `test=identity`, no exclusions beyond those the base-against-base comparison declares.
- Reported, not gated: the class averages of C0 and the adopted arm for the maintainer to inspect (the three sets,
  seed 1).

## 7. Common rules and stop conditions

- No floor relaxed; a file outside the table gets a row before the edit; the thread policy.
- Stops: the off-state identity breaks; the pose student changes; the alignment would need a change in
  `src/main/pftc` or `src/main/class`; C fails its wiring check after one fix; a symmetry the amplitude path would
  give up; decision 7 (c) not met (the table goes to the maintainer).
- Gates no ruling may waive: the off-state identity, the pose student's identity, thread identity of training.

## 8. Phases

- **Phase 0, the baseline.** Pin the base; Release and Debug builds, fast gate. C0 and R0 (R0 needs the validation
  to skip centring: a patch on an export for this phase, the same switch Phase 1 builds properly) at seeds 1 to 8 on
  the three sets; every read of decision 6 that exists without the alignment; decision 8 checked.
  Exit: the tables frozen in Progress.
- **Phase 1, alignment after prediction.** In the tree: `predict_classes` aligns every predicted particle (S1 and
  SK) and writes the committed class, angle, shift and route; `validate_class_student` drops `center_ptcls`, zeroes
  the unseen particles' shifts, measures both alignments and adds decision 6's reads; its report says what the
  shortlist is and which particles its class averages use, and names the alignment. Testers and plants.
  Exit: section 5's Phase 1 gates.
- **Phase 2, the alignment's read.** RS1 and RSK at seeds 1 to 8 on the three sets. Exit: decision 7 (a) applied;
  the alignment SK or S1 kept in the tree (the other removed).
- **Phase 3, the amplitude channel.** B, B', C as patches on an export of Phase 2's tree; C's wiring check; B and B'
  at seeds 1 to 8 on the three sets. Exit: decision 7 (b) applied.
- **Phase 4, adoption.** Decision 7 (c): the route without centring kept (B in the tree if adopted, its model
  version raised), `doc/algorithms/class_student.md` and the test policy updated; or the stop with the table.
  Exit: section 5's gates, the report.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90` (and tester) | prediction with alignment S1/SK | 1, 2, 4 |
| `src/main/strategies/parallelization/simple_class_student_align.f90` (and tester) | the alignment after prediction, a class of its own over the classic searches' shift map (row added in Phase 1, ruling R3) | 1, 2, 4 |
| `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` | the alignment tester's suite in `unit_numerics` and the suite help (row added in Phase 1, ruling R3) | 1, 2, 4 |
| `src/main/commanders/simple/simple_commanders_class_student.f90` | validation without centring, the reads, the report | 1, 2, 4 |
| `src/main/learn/simple_class_student.f90`, `simple_polar_student.f90`, `simple_polar_feats.f90`, `simple_hhcnn.f90` (and testers) | the amplitude path (if adopted) | 4 |
| `src/defs/simple_defs_hhcnn.f90` | named constants of the alignment and the amplitude path | 1, 2, 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` | the workflow test's reads and labels | 1, 2, 4 |
| `doc/algorithms/class_student.md`, `doc/policies/test_environment_policy.md` | the route and the test (Phase 2 added for the losing alignment's text, ruling R3) | 1, 2, 4 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/polar_student_future_directions.md` | record | 0, 1, 2, 3, 4 |
| `doc/implementation_notes/planned/pose_student_amplitude_channel.md` (deleted), `doc/implementation_notes/completed/class_student_vs_reference_matching.md`, `doc/implementation_notes/planned/selection_student.md`, `doc/implementation_notes/planned/student_inplane_posterior.md`, `doc/implementation_notes/planned/student_inplane_posterior_report.md` | the amplitude channel idea removed from the docs (rows added in Phase 4, ruling R4) | 4 |

## Progress

### Phase 0, the baseline (2026-11-16, the Dell, done)

**The base and the machine.** HEAD `1637eeb1ca1d26e6f5783a9d0e809841c86258b8` with the launch's uncommitted
documentation (this plan, untracked, and its row in `doc/implementation_notes/README.md`); no source file differs
from HEAD, so the source the run builds is `git archive 1637eeb1c`, exported to `scratch/keep/src_base`. The Dell:
Oracle Linux 8.10, Intel Core Ultra 9 285 (24 cores, no hyperthreading), 62 GB, GCC 15.2 of gcc-toolset-15, CMake
3.26. Run root `/data2/agent_runs/class_raw` (`R` below; paths under `R/scratch/keep/` unless absolute). The
maintainer's launch rulings apply: the run builds and runs on the Dell without following master, never uses the
maintainer's checkout or its build, and records the machine's load with every timing (ruling R1). Every script
points `SIMPLE_PATH` and `PATH` at the tree it means (`tools/env.sh <tree>`), since the login shell points them at
the maintainer's build.

**Builds.** Configured as `cmake -S src_base -B <tree> -DCMAKE_BUILD_TYPE=<Release|debug> -DBUILD_TESTS=ON
-DNICE=OFF -DUSE_ARCHOPT=ON -DFORTLS=OFF -DPython3_EXECUTABLE=/usr/bin/python3.12`, `make -j20`, `make install`
(`tools/build_base.sh`): `build_base` (Release, make 86 s) and `build_base_debug` (Debug, 42 s). The fast gate
(`scripts/run_fast_gate.sh`): 15 of 15 on both (9.7 s and 38.0 s; `R/scratch/logs/build_base{,_debug}.gate.log`).

**Off-state identity, the base against itself.** `simulated_workflow` (the end-to-end simulated test: solve2D's
refine2D stages, then solve3D's refine3D stages in shared memory and over two parts), both suites, at one thread
and `SIMPLE_SEED=7`, twice on `build_base` (`identity/base1`, `identity/base2`, side by side; 1jxy 14:13 and 13:56,
6vxx 14:02 and 12:16 minutes); `test=identity` of `base1` against `base2`: 1jxy 223 files the same and 0 differing,
6vxx 186 and 0 (`identity/cmp_base1_base2/<suite>/identity_report.txt`). **The exclusion lists are empty for both
suites.** 1jxy's shared-memory solve3D misses its pose gate in both runs (pair pose error 19.8° against 15°), the
known miss at this seed. `base1` stays for Phases 1 and 4; `base2`'s outputs were deleted, its logs kept.

**Decision 4.** `/data/10335-streptav/tst/1_import_particles/tst.simple` holds 11,402 particles, all selected, all
at zero shift: the maintainer's 11,402 (his validation's copy of the project is the same 11,402 rows). It is the
streptavidin set; `3_selection` (6,821 selected, with a solve2D run's shifts) is not used. bgal: 5,513 particles,
trpv1: 35,645, both all selected at zero shift.

**The arms.** One `validate_class_student` per set and seed (`tools/validate.sh`, `tools/chain_p0.sh`):
`SIMPLE_SEED=<seed> simple_exec prg=validate_class_student projfile=<copy of the import project> nsample=<n>
ncls=<k> mskdiam=<m> nthr=24`, with streptavidin `5000 30 85`, bgal `2500 50 180`, trpv1 `2500 50 175`, seeds 1 to
8. The validation draws `nsample` training particles, runs solve2D on them (the teacher), trains the student on the
teacher's class averages and posterior, labels the other particles (the unseen: 6,402 on streptavidin, 3,013 on
bgal, 33,145 on trpv1) and compares them with the reference.
- **C0**, the base as it is (`build_base`): every particle centred by `center_ptcls` first, then the student's
  prediction, then today's registration (a local gradient refinement of the shift, within 25 Å, at the predicted
  class and angle).
- **R0**, the base with one patch (`variants/R0/R0.patch`, 53 lines, on a copy of the export in
  `variants/R0/src`, built Release in `variants/R0/build`, fast gate 15 of 15): the validation skips
  `center_ptcls` and sets the unseen particles' shifts to zero (they already are, in all three import projects);
  everything else is the base's. This is the switch Phase 1 builds properly.
- The two arms ran back to back per set and seed, their order alternated by seed parity, nothing else heavy on the
  machine; the one-minute load at the start of a run was 10 to 22 (the runs' own 24 threads; 2.9 for the first).
  At the same seed the arms draw the same training particles (the same `teacher_selection_hash` in the model
  header), but solve2D at 24 threads does not reproduce its class averages, so the pairs share a seed and a draw,
  not a teacher.
- One run crashed: C0 bgal seed 1 in the chain, a segmentation fault two seconds in, entering solve2D's first
  refine2D stage (the backtrace: `sp_project`'s `write_segment_inside` called from the refine2D strategy's
  in-memory iteration, inside libc). The same tree, set and seed passed before (a smoke run) and after (the re-run
  the tables use), so it is not deterministic: a base fault outside this plan's files, left for the maintainer
  (log `validate/C0/bgal_s1_crash1/validate.log`). A Debug build (bounds-checked) of the same
  case passed as well (`debug_crash/C0_bgal_s1`, 223 s): no out-of-bounds access caught.

**What each read is** (`tools/reads.py`, from the validation's own outputs; every read is over the unseen
particles). The reference is the physical posterior's best class: one full search of the dense 2D probability
table over every class of the teacher, every in-plane rotation and the shift disk (three widths of a 10 Å prior)
at the full box, from the shifts the particles carry into it: in C0 the centred ones (the base's yardstick), in R0
the raw picks. The student's candidate list is the classes holding 99% of its posterior, at most 16; a particle
whose confidence reaches the model's calibrated threshold and whose list fits the cap is given a shortlist (the
only classes a seeded search tests). The mode is the student's most likely class. The shift and angle reads
compare the student's registered shift and read-out angle with the reference's (same class: where the mode is the
reference's class; all: every particle). The picking error is the length of the reference's shift (the raw picks
are at zero). The class averages are made from every unseen particle at the student's class, angle and shift (the
student's) or at the reference's (the reference's), and are compared with the teacher's on the classes holding at
least ten particles in both assignments: the FRC 0.5 crossing (Fourier ring correlation) interpolated between
shells, and the mean FRC over the shells to 15 Å. Costs are wall time per particle at 24 threads: centring over
every particle, prediction (features, forward pass and registration, one timer in the base) and the full search
over the unseen. `reads.py` reproduces the validation's own shortlist, hold, agreement, shift and FRC values. Its
interpolated crossing starts its walk where the validation's shell read does: a fix made during the run, since the
copy from the last run returned the box size when the first shell lay below 0.5; the five rows written before the
fix were recomputed from the kept stacks.

**The table** (means over the eight paired seeds; R0 - C0 with its paired standard error, * beyond two; rows in
`R/scratch/keep/reads.tsv`, this table in `progress_table_p0.md`, the per-set tables with more reads in
`tables_p0_C0_R0.md`):

| read | strep C0 | strep R0 | R0 - C0 (SE) | bgal C0 | bgal R0 | R0 - C0 (SE) | trpv1 C0 | trpv1 R0 | R0 - C0 (SE) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| particles given a shortlist (%) | 3.9 | 0.7 | -3.2 (0.7) * | 85.6 | 54.7 | -31.0 (2.6) * | 91.8 | 64.2 | -27.6 (2.0) * |
| shortlists holding the reference's class (%) | 99.4 | 94.2 | -5.2 (1.8) * | 92.8 | 94.5 | +1.7 (0.7) * | 97.5 | 98.1 | +0.6 (0.3) |
| shortlist length (median of the run's medians) | 14.4 | 14.4 | +0.0 (0.4) | 3.9 | 7.8 | +3.9 (0.3) * | 4.2 | 9.2 | +5.0 (0.2) * |
| candidate list holds the reference's class, all particles (%) | 94.4 | 92.5 | -1.9 (0.5) * | 91.3 | 89.5 | -1.8 (0.7) * | 97.1 | 96.4 | -0.7 (0.3) * |
|   same, the nearer half of the picking error (%) | 96.2 | 96.5 | +0.3 (0.3) | 94.3 | 97.7 | +3.4 (0.5) * | 98.2 | 99.4 | +1.3 (0.2) * |
|   same, the farther half (%) | 92.7 | 88.6 | -4.2 (0.7) * | 88.4 | 81.3 | -7.1 (1.0) * | 96.1 | 93.4 | -2.7 (0.5) * |
| candidate list length (mean) | 15.91 | 15.98 | +0.08 (0.02) * | 6.75 | 11.78 | +5.03 (0.24) * | 6.26 | 11.76 | +5.50 (0.27) * |
| the mode is the reference's class (%) | 28.1 | 22.6 | -5.6 (0.8) * | 48.9 | 33.2 | -15.7 (0.9) * | 54.5 | 35.3 | -19.3 (1.0) * |
|   same, predicted at the reference's shifts (%) | 39.1 | 37.5 | -1.5 (1.4) | 65.7 | 53.3 | -12.4 (1.0) * | 72.1 | 59.9 | -12.1 (0.8) * |
| shift to the reference's, same class, median (Å) | 1.39 | 2.08 | +0.69 (0.07) * | 0.17 | 0.29 | +0.13 (0.01) * | 0.25 | 0.59 | +0.34 (0.01) * |
| shift to the reference's, same class, 90th pct (Å) | 6.33 | 11.66 | +5.33 (0.25) * | 1.94 | 7.96 | +6.02 (0.58) * | 1.65 | 7.43 | +5.78 (0.24) * |
| shift to the reference's, all, median (Å) | 3.08 | 6.56 | +3.48 (0.19) * | 1.41 | 5.43 | +4.02 (0.20) * | 0.99 | 3.79 | +2.79 (0.08) * |
| shift to the reference's, all, 90th pct (Å) | 12.92 | 18.33 | +5.41 (0.42) * | 12.09 | 22.69 | +10.60 (0.25) * | 6.79 | 14.07 | +7.28 (0.15) * |
| in-plane angle to the reference's, same class, median (deg) | 2.95 | 6.04 | +3.10 (0.55) * | 1.08 | 1.81 | +0.73 (0.04) * | 1.14 | 2.29 | +1.15 (0.03) * |
| picking error: the reference's shift, median (Å) | 10.00 | 7.26 | -2.74 (0.24) * | 9.51 | 8.78 | -0.73 (0.13) * | 7.64 | 7.25 | -0.39 (0.05) * |
| picking error: the reference's shift, 90th pct (Å) | 19.66 | 17.04 | -2.62 (0.19) * | 20.82 | 21.90 | +1.08 (0.24) * | 14.93 | 15.13 | +0.20 (0.05) * |
| the shift entering prediction to the reference's, median (Å) | 3.32 | 7.26 | +3.94 (0.21) * | 5.40 | 8.78 | +3.39 (0.23) * | 2.86 | 7.25 | +4.39 (0.12) * |
| classes compared (FRC) | 27.9 | 27.8 | -0.1 (0.5) | 36.0 | 38.0 | +2.0 (1.2) | 26.4 | 31.9 | +5.5 (1.0) * |
| student's class averages, FRC 0.5 interpolated (Å) | 13.20 | 13.50 | +0.31 (0.10) * | 16.39 | 19.65 | +3.26 (0.68) * | 23.08 | 23.68 | +0.60 (0.99) |
| reference's class averages, FRC 0.5 interpolated (Å) | 12.56 | 12.66 | +0.11 (0.11) | 13.76 | 14.67 | +0.91 (0.38) * | 18.47 | 17.12 | -1.35 (0.57) * |
| student's class averages, mean FRC to 15 Å | 0.871 | 0.843 | -0.028 (0.011) * | 0.777 | 0.771 | -0.006 (0.014) | 0.703 | 0.727 | +0.024 (0.011) * |
| reference's class averages, mean FRC to 15 Å | 0.914 | 0.912 | -0.002 (0.009) | 0.825 | 0.811 | -0.014 (0.010) | 0.786 | 0.794 | +0.008 (0.005) |
| cost: centring (ms per particle) | 8.773 | 0.000 | -8.773 (0.060) * | 2.636 | 0.000 | -2.636 (0.036) * | 1.588 | 0.000 | -1.588 (0.046) * |
| cost: prediction with registration (ms) | 0.150 | 0.163 | +0.013 (0.003) * | 0.617 | 0.649 | +0.032 (0.003) * | 0.395 | 0.402 | +0.007 (0.005) |
| cost: the student's route, sum (ms) | 8.923 | 0.163 | -8.760 (0.058) * | 3.253 | 0.649 | -2.604 (0.035) * | 1.983 | 0.402 | -1.581 (0.044) * |
| cost: the reference's full search (ms) | 0.421 | 0.440 | +0.019 (0.004) * | 1.123 | 1.280 | +0.157 (0.004) * | 0.862 | 0.915 | +0.053 (0.005) * |
| validation wall time (s) | 164.9 | 65.5 | -99.4 (1.3) * | 68.8 | 54.2 | -14.6 (0.4) * | 180.7 | 124.0 | -56.7 (2.0) * |

**What it shows.**
- **The candidate list survives the raw picks.** It holds the reference's class for 92.5% (streptavidin), 89.5%
  (bgal) and 96.4% (trpv1) of the raw picks, 0.7 to 1.9 points below C0. The loss is in the farther half of the
  picking error (3 to 7 points), while the nearer half holds it as well or better. The lists grow from about 6 to
  12 classes on average on bgal and trpv1 (streptavidin's are at the cap of 16 of 30 in both arms).
- **The mode does not.** It agrees with the reference for 22.6, 33.2 and 35.3% of the raw picks (C0 28.1, 48.9
  and 54.5%): a search over the list, not the mode alone, is what the raw picks need.
- **The shortlist shrinks and holds.** 55% (bgal) and 64% (trpv1) of the raw picks are shortlisted against 86 and
  92% centred, with lists twice as long, and the shortlists hold the reference's class as often (94.5 and 98.1%).
  On streptavidin almost no particle is shortlisted in either arm (0.7 and 3.9%).
- **Today's registration does not centre raw picks.** The registered shift misses the reference's by 3.8 to 6.6 Å
  at the median and 14 to 23 Å at the 90th percentile, about the picking error itself (7.3 to 8.8 Å median, 15 to
  22 Å at the 90th percentile); centred, it misses by 1.0 to 3.1 Å.
- **The class averages hold or nearly hold.** The mean FRC to 15 Å of the student's averages is lower on
  streptavidin (0.843 against 0.871, beyond two standard errors), equal on bgal and higher on trpv1 (0.727
  against 0.703); the interpolated FRC 0.5 crossing is worse on streptavidin (+0.3 Å) and bgal (+3.3 Å). The
  reference's own averages move less (mean FRC within 0.014, the crossing within 1.4 Å).
- **Cost.** Without centring the student's route costs 0.16, 0.65 and 0.40 ms per particle, below the reference's
  full search (0.44, 1.28, 0.92 ms); with it, 8.9, 3.3 and 2.0 ms. A validation takes 65, 54 and 124 s instead of
  165, 69 and 181 s.

**Decision 8, checked.** R0's candidate list holds the reference's class for at least 88% of the particles in
every run on every set, far above half: the early stop does not apply, and RSK (the alignment over the whole list)
runs in Phase 2.

**Exit.** C0 and R0 at eight seeds on the three sets, the reads frozen above (`R/scratch/keep/reads.tsv`, 48
rows; the runs under `validate/<arm>/<set>_s<seed>/`). Files changed in this phase: this plan (the notes under
decisions 4 and 8, the status line, this entry). Kept for later phases: `build_base`, `build_base_debug`,
`variants/R0` (source, patch, build), `identity/base1`, the tools, seed 1's class-average stacks of both arms
(`validate/<arm>/<set>_s1/keep_cavgs/`), and every run's report, details, log, model header and prediction tables.

### Phase 1, the alignment after prediction (2026-11-16, the Dell, done)

**What was built.** Given the teacher's class averages, `predict_classes` now aligns every particle after the
prediction, at its picked position: every in-plane rotation of the polar grid and shifts over the disk of the
validation's reference pass (radius three widths of the 10 Å prior, 30 Å; under that prior; at the student's band
and the physical temperature). The local shift registration of the base (a gradient refinement within 25 Å at the
predicted class and angle) is gone from the class student; the pose student keeps its own. The search is reused
unchanged: the classic 2D and 3D searches' polar shift map object (`shiftmap_srch` in `strategies/search`, which
calls the calculator's coarse shift map and the gradient shift search's refinement from a map cell, both in
`src/main/pftc`, untouched). Per class: the score at the stored shift at its best rotation, the coarse map over
the disk, the refinement from the class's best cell; a class keeps the refined pose when it beats its stored-shift
score; the best class by score is kept, as the classic 2D search commits. It lives in a small class of its own,
`class_student_align` (`src/main/strategies/parallelization/simple_class_student_align.f90`, one object per thread).
The two variants run on the same prediction in separate parallel loops, each timed: **S1** against the most likely
class, **SK** against every class of the candidate list (at most 16). In this phase SK's class, angle and shift are
committed to the project; S1's are written beside them in the prediction table (`class_s1 e3_s1 x_s1 y_s1`), and the
wall times of the prediction and of each alignment go to `class_student_times.txt` (a shared-memory run's). Before
the alignment the particles are marked not yet centred (`shsig` −1), after it their centring state records the move
(`record_shift_state`), as a 2D search or `center_ptcls` leaves it. `predict_classes` given class averages sets the
alignment's prior width (`shift_sigma` 10 Å, `CLASS_DENSE_SHIFT_SIGMA`), the window (`trs` 5 px, unused by the
disk) and the physical temperature (`prob_beta` 1); no new command-line key.

`validate_class_student` no longer runs `center_ptcls`: the unseen particles enter at their picked positions (shift
zero, not yet centred), the same project goes to the prediction and to the reference pass. It reads both alignments
from one prediction, makes the class averages of each (`cavgs/cavgs_student.mrc` for SK, `cavgs_s1/` for S1) and
reports decision 6's reads: in the report the shortlist (what it is, its share, its hold on the reference's class,
its length), each alignment's agreement with the reference's class, each alignment's class averages against the
teacher's (FRC 0.5 interpolated between shells, beside the reference's), and the route's cost per particle
(prediction plus one alignment) against the full search; in the details the network's reads at the picked positions
and at the reference's shifts, per alignment the shift's distance to the reference's (same class and all, median
and 90th percentile, Å) and the in-plane angle's, the picking error, the class averages' FRC 0.5 (whole shells and
interpolated) and mean FRC to 15 Å on one common class set, and the costs of the prediction and of each alignment.

**Decisions taken under delegation.**
- The alignment reuses `shiftmap_srch` (the classic searches' object) rather than the dense table's inline fill code:
  the same calculator routines (coarse map, mode, cell refinement), wrapped once already for per-particle use. The
  commit rule is the classic 2D search's (the refined pose when it beats the stored shift's score; the best class by
  score), so S1 and SK differ only in the list. The map's mode is taken under the prior at the particle's
  temperature, as the reference pass treats a particle not yet centred.
- SK is committed in this phase (the candidate list includes the mode, so SK is S1 or better by score); Phase 2
  keeps one alignment.
- A new tester suite `class_student_alignment` in `unit_numerics`, which needed rows for the suite table and the
  test UI's suite help (ruling R3, rows added to the file table with this phase).
- The report's reference-line wording changed to "The reference is ..." so that the workflow test's label for the
  reference pass's run time matches the timing row, not the explanation.
- The recorded angle is wrapped into [0°, 360°) (`modulo`), the student's own read-out convention: `360 − rot`
  gave 360.0 at the first rotation, which the workflow test's record gate refused on its first run.

**Testers and plants.** `simple_class_student_align_tester.f90` (suite `class_student_alignment`, 172 checks): white
class averages on the band, each particle one of them rotated by a grid rotation and shifted by 5 px (beyond the 2 px
window, inside the 6 px disk): S1 recovers the planted class, rotation and shift (within 0.25 px; the increment is
the planted shift itself); SK keeps the planted class wherever the list holds it, and another class alone aligns
worse; a team of four aligns as one thread, bit for bit (class, rotation, score, shift, width). Every check was
broken on purpose first, on a copy of the tree (`R/scratch/keep/tools/plants_p1.py`, patches and log in
`R/scratch/keep/plants/`): the refinement ignored (135 failures), the first class kept (81), the window instead of
the disk (124), a thread-dependent score (6, the team checks); the clean and the restored copies pass. On the way the
fixture showed that a test's sigma2 must outlive the fixture routine (the calculator keeps a pointer to it).

**Gates.**
- Builds of the run tree configured as Phase 0's (`repo/build` Release, `repo/build_debug` Debug). Fast gate 15 of
  15 on both on the final source, `unit_numerics` and `unit_pftc_align2D3D` among them (Release 6.6 s and 9.7 s;
  `R/scratch/logs/p1_gate_release_final.log`, `p1_gate_debug_final.log`).
- CTest `class_student`: all ten gates pass (record, parts against shared memory byte for byte, the report's reads,
  agreement above chance) after the angle fix; its reported reads in `R/scratch/keep/ctest_p1/`. CTest
  `pose_student`: pass (36 s).
- Off-state identity: `simulated_workflow` at one thread, `SIMPLE_SEED=7`, on the run tree (`identity/tree`) against
  `base1`: 1jxy 223 files the same and 0 differing, 6vxx 186 and 0, no exclusions
  (`identity/cmp_base1_tree/<suite>/identity_report.txt`); 1jxy's known pose-gate miss with the base's value. The
  one source change after these runs and the guard below, the angle's wrap, lies in `predict_classes` alone, which
  neither runs.
- The pose student byte-identical: `learn_pose_student` on bgal at seed 1 and 24 threads (the student_speed set's
  teacher map) on `build_base` and on the run tree: the same header (paths and times excepted) and weights (md5
  `8f0d176e…`); the class student's training too (student_inplane's frozen bgal teacher: same header, weights
  `3f4517de…`) (`tools/guard_p1.sh`, `R/scratch/keep/guard_p1/`).
- Thread identity of the new loops, by hand: `predict_classes` with the alignment on bgal's 3,013 unseen particles at
  24 threads twice and at one thread, seeded: tables and candidates byte-identical. Unseeded, two 24-thread runs
  differed in one particle's shift (third decimal): the shift map plans its transforms with `FFTW_MEASURE` unless a
  global seed is set (`simple_polarft_shiftmap.f90`), existing behaviour of the reused code; every validation and
  test runs seeded (`simple_test_exec` fixes `SIMPLE_SEED`).

**The validation in one run, seed 1** (24 threads, the machine otherwise idle, one-minute load 4 to 17 from the run
itself; `validate/P1/<set>_s1`, rows in `R/scratch/keep/reads2.tsv` by `tools/reads2.py`, which reproduces the
validation's own values; one seed: a check that the route works, not a verdict):

| read (seed 1) | streptavidin | bgal | trpv1 |
| --- | --- | --- | --- |
| candidate list holds the reference's class (%) | 89.9 | 88.9 | 96.4 |
| the most likely class is the reference's (%) | 22.5 | 32.9 | 32.7 |
| S1's aligned class is the reference's (%) | 22.5 | 32.9 | 32.7 |
| SK's aligned class is the reference's (%) | 28.8 | 54.6 | 69.7 |
| SK's shift to the reference's, all, median / 90th pct (Å) | 5.08 / 16.9 | 0.00 / 14.8 | 0.00 / 6.4 |
| S1's shift to the reference's, all, median / 90th pct (Å) | 5.91 / 18.0 | 3.25 / 21.5 | 2.67 / 12.3 |
| class averages, mean FRC to 15 Å: S1 / SK / reference | 0.850 / 0.905 / 0.894 | 0.822 / 0.840 / 0.847 | 0.771 / 0.776 / 0.819 |
| class averages, FRC 0.5 interpolated (Å): S1 / SK / reference | 13.2 / 12.3 / 12.5 | 18.5 / 14.2 / 14.3 | 20.5 / 14.5 / 14.6 |
| cost per particle (ms): prediction / S1 / SK | 0.154 / 0.021 / 0.354 | 0.665 / 0.041 / 0.532 | 0.434 / 0.035 / 0.466 |
| route with S1 / with SK / full search (ms) | 0.18 / 0.51 / 0.42 | 0.71 / 1.20 / 1.30 | 0.47 / 0.90 / 0.91 |
| validation wall time (s) | 64 | 56 | 157 |

Where the aligned class is the reference's, the aligned pose is the reference's to the rounding (shift within
0.006 Å, the same rotation): the alignment and the reference pass run the same map and refinement from the same
start. The prediction's cost includes reading and preparing the particles.

**Files changed in this phase:** `src/main/strategies/parallelization/simple_strategy_class_student.f90` (the
prediction with both alignments), `simple_class_student_align.f90` and `simple_class_student_align_tester.f90` (new,
103 and 270 lines), `src/main/commanders/simple/simple_commanders_class_student.f90` (`predict_classes`' settings,
the validation), `src/defs/simple_defs_hhcnn.f90` (the times file and S1's class-average stack names; the prior's
comment), `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` (the reads and labels),
`src/main/commanders/test/simple_commanders_test_class.f90` and `src/main/ui/simple_test/simple_test_ui_class.f90`
(the suite), `doc/algorithms/class_student.md`, `doc/policies/test_environment_policy.md`,
`doc/implementation_notes/planned/polar_student_future_directions.md` (section 6, the candidate list as a
neighbourhood for the pose student's post-prediction refinement), this plan (file-table rows, this entry).

**Exit.** Section 5's Phase 1 gates: the testers of the alignment with their plants, the fast gate,
`unit_numerics`, `unit_pftc_align2D3D` on both trees, the `class_student` workflow test, the off-state identity;
and the pose student's identity. Kept for Phase 2: `repo/build`, `tools/validate2.sh`, `tools/reads2.py`,
`reads2.tsv`.

### Phase 2, the alignment's read (2026-11-16, the Dell, done)

**What was run.** RS1 and RSK, the student on raw picks with the alignment against its most likely class (S1) or
against every class of its candidate list (SK), come from one `validate_class_student` run each: Phase 1's tree runs
both alignments on the same prediction. Seeds 1 to 8 on the three sets with Phase 0's command lines (seed 1 is
Phase 1's run of the same source; seeds 2 to 8 `tools/chain_p2.sh` with `tools/validate2.sh`, arm label `P1`,
24 threads, nothing else heavy on the machine, one-minute load at the start of a run 2 to 23, the runs' own
threads; 24 runs, no failure). Rows in `R/scratch/keep/reads2.tsv` (`tools/reads2.py`, which cross-checks the
validation's own values); the paired tables with every read in `R/scratch/keep/tables_p2_S1_SK.md`
(`tools/p2_table.py`), this one in `progress_table_p2.md`.

**The table** (means over the eight seeds; SK − S1 paired within each run, its standard error, * beyond two):

| read | strep S1 | strep SK | SK - S1 (SE) | bgal S1 | bgal SK | SK - S1 (SE) | trpv1 S1 | trpv1 SK | SK - S1 (SE) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned class is the reference's (%) | 22.7 | 30.0 | +7.4 (0.4) * | 34.0 | 53.2 | +19.3 (1.0) * | 34.1 | 72.1 | +37.9 (0.5) * |
| shift to the reference's, all, median (Å) | 5.45 | 4.98 | -0.47 (0.15) * | 2.94 | 0.00 | -2.94 (0.12) * | 2.36 | 0.00 | -2.36 (0.12) * |
| shift to the reference's, all, 90th pct (Å) | 16.83 | 15.91 | -0.92 (0.19) * | 21.86 | 18.08 | -3.77 (0.94) * | 12.64 | 6.99 | -5.65 (0.19) * |
| class averages, FRC 0.5 interpolated (Å) | 13.08 | 12.65 | -0.43 (0.13) * | 16.12 | 13.81 | -2.31 (0.74) * | 21.84 | 17.57 | -4.26 (1.14) * |
| class averages, mean FRC to 15 Å | 0.874 | 0.913 | +0.039 (0.007) * | 0.808 | 0.832 | +0.024 (0.006) * | 0.764 | 0.780 | +0.016 (0.011) |
| alignment's cost (ms per particle) | 0.024 | 0.388 | +0.364 (0.007) * | 0.040 | 0.534 | +0.494 (0.010) * | 0.036 | 0.448 | +0.412 (0.006) * |
| route's cost, prediction and alignment (ms) | 0.179 | 0.543 | +0.364 (0.007) * | 0.693 | 1.186 | +0.494 (0.010) * | 0.439 | 0.850 | +0.412 (0.006) * |

| beside them (means of the same runs) | strep | bgal | trpv1 |
| --- | --- | --- | --- |
| candidate list holds the reference's class (%) | 92.2 | 90.1 | 96.3 |
| SK keeps the most likely class (%) | 12.8 | 29.1 | 31.3 |
| reference's class averages, FRC 0.5 interpolated (Å) | 12.73 | 14.16 | 17.30 |
| reference's class averages, mean FRC to 15 Å | 0.913 | 0.832 | 0.806 |
| full search (ms per particle) | 0.445 | 1.276 | 0.902 |
| validation wall time (s) | 67.5 | 57.7 | 151.4 |
| C0 (base, centred): the committed class is the reference's (%) | 28.1 | 48.9 | 54.5 |
| C0: class averages, mean FRC to 15 Å | 0.871 | 0.777 | 0.703 |
| C0: route's cost (ms) | 8.923 | 3.253 | 1.983 |

Read with Phase 0's columns: "aligned class" is the class committed to the project; for C0 and R0 it was the most
likely class. Where the aligned class is the reference's, the aligned pose is the reference's to the rounding (the
same map and refinement from the same start), so the shift's median is 0 Å wherever more than half the particles
agree. The validation's per-particle costs include reading and preparing the particles.

**Decision 7 (a): SK.** Its aligned class agrees with the reference's more often than S1's on every set, far beyond
twice the paired standard error (+7.4, +19.3 and +37.9 points); the cheaper S1 is chosen only within that margin.
SK also brings the shifts closer to the reference's and gives class averages nearer the teacher's (FRC 0.5 crossing
0.4 to 4.3 Å finer; mean FRC to 15 Å higher on streptavidin and bgal, equal within its error on trpv1), at 0.36 to
0.49 ms more per particle. On raw picks SK's agreement now exceeds the centred base's (C0) on bgal and trpv1 (53.2
against 48.9, 72.1 against 54.5) and nearly matches it on streptavidin (30.0 against 28.1). For decision 7 (c), to
be judged in Phase 4: SK's route costs 0.54 ms on streptavidin against the full search's 0.45 ms (above it), 1.19
against 1.28 on bgal and 0.85 against 0.90 on trpv1 (below).

**The losing alignment removed.** S1 is gone from the tree: `predict_classes` runs one alignment against the
candidate list (`align_one`), its prediction table is back to the base's columns, its times file holds the
prediction and the alignment; the validation reads one alignment (report: the aligned class's agreement, the class
averages' FRC 0.5 crossing beside the reference's, the route's cost; details: the shifts, the angle, the picking
error, the class averages' FRC on the student's and the reference's common classes, the costs); the S1 class-average
stack's name is gone from `simple_defs_hhcnn.f90`; the workflow test reads 18 values; the tester's checks keep
their content (a list of one class recovers the planted pose; a list of every class keeps the planted class) with
the variants' names dropped from their wording; the algorithm page and the test policy describe the one route (the
row of those two documents gains this phase, ruling R3). The alignment module and its checks are unchanged, so
Phase 1's plants still apply. Decision under delegation: the details' class-average block compares the student's
and the reference's averages on the classes both assignments populate (Phase 1 used the three assignments').

**Phase 1's gates again**, on the source without S1:
- Fast gate 15 of 15 on both trees, `unit_numerics` and `unit_pftc_align2D3D` among them (`R/scratch/logs/
  p2_gate_release_final.log`, `p2_gate_debug_final.log`).
- CTest `class_student`: all ten gates pass (`R/scratch/keep/ctest_p2/class_student_metrics.tsv`); `pose_student`:
  pass.
- Off-state identity, `simulated_workflow` at one thread against `base1`: 1jxy 223 files the same and 0 differing,
  6vxx 186 and 0 (`identity/cmp_base1_tree_p2/`).
- The pose student and the class student's training byte-identical to the base on bgal at seed 1 (headers and
  weights, `tools/guard_p2.sh`, `R/scratch/keep/guard_p2/`).
- Seeded thread identity of the prediction with the alignment, by hand: bgal's 3,013 unseen particles at 24 threads
  twice and at one thread give the same table and candidates byte for byte (`R/scratch/p2/stress/`).
- One validation of the final tree (bgal, seed 1, `validate/P2/bgal_s1`, `tools/validate3.sh` with `tools/reads3.py`
  for the later phases): the aligned class agrees for 52.5%, the route costs 1.03 ms against the full search's 1.14 ms.

**Files changed in this phase:** `src/main/strategies/parallelization/simple_strategy_class_student.f90`,
`src/main/commanders/simple/simple_commanders_class_student.f90`, `src/defs/simple_defs_hhcnn.f90`,
`src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90`,
`src/main/strategies/parallelization/simple_class_student_align_tester.f90`, `doc/algorithms/class_student.md`,
`doc/policies/test_environment_policy.md`, `doc/implementation_notes/planned/polar_student_future_directions.md`
(the verdict noted), this plan (the documents' row, this entry).

**Exit.** Decision 7 (a) applied: SK kept in the tree, S1 removed; Phase 1's gates pass on the result.

### Phase 3, the amplitude channel (2026-11-16, the Dell, done)

B, B' and C (decision 3) were built as patches on an export of Phase 2's tree, outside the tree, with the network
tester's gradient, symmetry and thread-identity checks extended to the amplitude path (132 checks pass; the pose
student's training byte-identical to the base's), and run beside the tree (arm P3, the raw-pick route with SK) at
seeds 1 to 8 on the three sets, 72 runs (every read in `R/scratch/keep/tables_p3_B_Bp.md`, the patches in
`R/scratch/keep/variants/`). C passed the wiring check as revised (decision 3): its reads did not move with the
shift, its angle could not resolve the half turn, and it recognised 0.484 of bgal's held-out classes. Decision 7
(b): B rejected, worse than the tree on reads of every set (streptavidin's aligned class −4.7 points and class
averages' mean FRC −0.012, bgal's class averages 0.65 Å coarser at FRC 0.5, trpv1's share shortlisted −10.0
points), its one large gain, bgal's aligned class (+8.0 points), shown by B' too (+6.0): capacity, not the
amplitudes; it cost 0.02 to 0.10 ms more per particle. The tree was unchanged. The amplitude channel was tested on
the class student 2026-11-16 and rejected; not to be revisited.

### Phase 4, the verdict on dropping centring and the adoption (2026-11-16, the Dell, done)

**Decision 7 (c), judged on the tables of Phases 0, 2 and 3.** The best raw-pick arm is the tree as it stands: the
student on raw picks with the alignment against its candidate list (SK, decision 7 (a)), without the amplitude path
(rejected by decision 7 (b)). It is compared with C0, the base with `center_ptcls` (Phase 0), seed by seed, in the
tree's own runs of Phase 3 and, as a replicate, the SK columns of Phase 2's runs (eight seeds each; differences with
their paired standard errors; "better" and "worse" beyond two of them; `R/scratch/keep/table_7c.md`):

| read | set | C0 (base, centred) | raw picks, SK (Phase 3) | difference (SE) | raw picks, SK (Phase 2) | difference (SE) |
| --- | --- | --- | --- | --- | --- | --- |
| committed class is the reference's (%) | strep | 28.1 | 30.7 | +2.6 (0.8) better | 30.0 | +1.9 (0.5) better |
|  | bgal | 48.9 | 53.3 | +4.4 (0.9) better | 53.2 | +4.3 (0.5) better |
|  | trpv1 | 54.5 | 71.9 | +17.4 (0.9) better | 72.1 | +17.5 (0.9) better |
| class averages, FRC 0.5 interpolated (Å) | strep | 13.20 | 12.45 | -0.75 (0.13) better | 12.65 | -0.55 (0.14) better |
|  | bgal | 16.39 | 13.46 | -2.93 (0.73) better | 13.81 | -2.58 (0.72) better |
|  | trpv1 | 23.08 | 18.28 | -4.80 (1.27) better | 17.57 | -5.51 (0.98) better |
| class averages, mean FRC to 15 Å | strep | 0.871 | 0.914 | +0.043 (0.012) better | 0.913 | +0.041 (0.009) better |
|  | bgal | 0.777 | 0.829 | +0.052 (0.010) better | 0.832 | +0.055 (0.013) better |
|  | trpv1 | 0.703 | 0.755 | +0.052 (0.021) better | 0.780 | +0.076 (0.017) better |
| cost of the route, prediction and alignment (ms per particle) | strep | 8.923 | 0.537 | -8.387 (0.054) better | 0.543 | -8.381 (0.054) better |
|  | bgal | 3.253 | 1.162 | -2.091 (0.028) better | 1.186 | -2.066 (0.030) better |
|  | trpv1 | 1.983 | 0.853 | -1.130 (0.045) better | 0.850 | -1.133 (0.054) better |
| cost of the full search, the same runs (ms per particle) | strep | 0.421 | 0.442 | +0.021 (0.009) | 0.445 | +0.025 (0.007) |
|  | bgal | 1.123 | 1.271 | +0.148 (0.008) | 1.276 | +0.152 (0.010) |
|  | trpv1 | 0.862 | 0.923 | +0.061 (0.021) | 0.902 | +0.040 (0.010) |

| route below the full search (paired within each run) | set | Phase 3: route - full (SE) | runs below | Phase 2: route - full (SE) | runs below |
| --- | --- | --- | --- | --- | --- |
| | strep | +0.095 (0.004) | 0 of 8 | +0.097 (0.003) | 0 of 8 |
| | bgal | -0.110 (0.014) | 8 of 8 | -0.089 (0.008) | 8 of 8 |
| | trpv1 | -0.070 (0.016) | 8 of 8 | -0.051 (0.008) | 8 of 8 |

- **Committed-class agreement and class averages: not worse than C0 on any set; better on all three.** The class a
  particle ends in agrees with the reference's 2.6, 4.4 and 17.4 points more often than C0's, and the class
  averages' FRC 0.5 crossing (Fourier ring correlation against the teacher's averages) is 0.75, 2.93 and 4.80 Å finer,
  their mean FRC to 15 Å higher by 0.043 to 0.052.
- **Cost: below the full search on bgal and trpv1, not on streptavidin.** The route (prediction and alignment)
  costs 0.537 ms per particle on streptavidin against the full search's 0.442 ms in the same runs: above it in every
  one of the eight runs, and in the eight of Phase 2 too. On bgal (1.162 against 1.271 ms) and trpv1 (0.853 against
  0.923 ms) it is below in every run. Against C0's route, which pays for the centring, it is 16, 2.8 and 2.3 times
  cheaper.
- **Verdict as written: not met** (cost on streptavidin). The rule asks for every set and ruling R2 keeps decision
  7's rules out of the run's reach, so the run stopped with the table for the maintainer (section 7).
- **Why streptavidin.** With 30 classes the candidate list is capped at 16 (`KMAX_CLASS`): SK aligns each particle
  against about half of all classes (15.98 on average, Phase 2), at 0.39 ms, while the full search over all 30 costs
  0.44 ms including its reading and preparation. On bgal and trpv1 (50 classes) the lists average 11.8 and 12.7.
  S1 would have been below (0.18 ms) but lost decision 7 (a) by 7 to 38 points of agreement.

**Ruling R4 (Hans, 2026-11-16), the adoption.** The cost condition of decision 7 (c) is waived for streptavidin, a
30-class set where the full search is already cheap; the raw-pick route is adopted as built (no centring, the
alignment after prediction against the candidate list, no amplitude path); C's revised wiring check is confirmed
(decision 3); the amplitude channel idea leaves the docs. A list cap relative to the class count stays a possible
follow-up, not done. Under it:
- **The amplitude channel removed from the docs.** `planned/pose_student_amplitude_channel.md` deleted and its
  README row dropped; the amplitude proposal and every pointer to that plan removed from
  `planned/polar_student_future_directions.md` (the moduli front end among section 2's alternatives, the planned
  path, section 6's note), `completed/class_student_vs_reference_matching.md` (section 7's paragraph after the NUViT
  table, the table kept), `planned/selection_student.md`, `planned/student_inplane_posterior.md` (its sentences
  about the amplitude plan and the shared Phase 0, fourteen places) and its report (two cells of its phase table;
  the `git status` listing it copied is kept verbatim). Phase 3's entry above is now a short paragraph; decision 3
  and the arms B, B', C stay as the record of what was run. The five files got a file-table row first.
- **The route's documents.** `doc/algorithms/class_student.md` already describes the route without centring, the
  alignment over the candidate list and the validation's reads (Phases 1 and 2); it gained the raw-pick reads under
  "What the reads found". The test policy's rows (`unit_numerics`, `class_student`) were right. The README row says
  adopted.
- **The gates of section 5 on the final tree** (source unchanged since Phase 2; the rebuild recompiled only the
  version string): the fast gate 15 of 15 on Release and Debug, with `unit_numerics` (the alignment tester's 172
  checks) and `unit_pftc_align2D3D` (`R/scratch/logs/p4r_gate_build*.log`); CTest `class_student` and
  `pose_student` pass (`p4_ctest_*.log`); the plants of Phase 1 stand (the alignment module is byte-identical to the
  copy they were planted in; the tester only renamed since); off-state identity, `simulated_workflow` 1jxy and 6vxx
  at one thread and `SIMPLE_SEED=7` on the final tree against the base: 223 and 186 files the same, 0 differing
  (`R/scratch/keep/identity/cmp_base1_tree_p4/`).
- **`validate_class_student` at seed 1 on the three sets with the final tree** (24 threads, alone on the machine,
  `tools/validate3.sh P4`, rows in `reads3.tsv`; C0 is the base's seed-1 run of Phase 0; one seed, so read beside
  the eight-seed table above; the 24-thread solve2D teacher is not the same twice, so these are not Phase 3's seed-1
  numbers again):

| read (seed 1) | streptavidin: route | C0 | bgal: route | C0 | trpv1: route | C0 |
| --- | --- | --- | --- | --- | --- | --- |
| committed class is the reference's (%) | 28.8 | 28.8 | 53.2 | 52.4 | 72.9 | 52.1 |
| candidate list holds the reference's class (%) | 93.4 | 94.8 | 92.0 | 93.6 | 97.3 | 97.9 |
| class averages, FRC 0.5 interpolated (Å); the reference's | 12.50; 13.05 | 13.16 | 13.65; 14.14 | 18.61 | 14.68; 21.13 | 24.23 |
| class averages, mean FRC to 15 Å | 0.897 | 0.884 | 0.824 | 0.778 | 0.788 | 0.696 |
| route's cost (ms per particle); the full search's | 0.484; 0.390 | 8.695 | 1.132; 1.235 | 3.372 | 0.792; 0.842 | 1.688 |
| validation wall time (s) | 56 | | 53 | | 130 | |

  The reports (`R/scratch/keep/validate/P4/<set>_s1/1_validate_class_student/`) name the alignment and mark the
  cost met on bgal and trpv1, not on streptavidin (waived). The shortlists hold the reference's class for 96.6% on
  bgal and 99.3% on trpv1, met; on streptavidin the student shortlists 9 of 6,402 particles (0.1%) and 8 of the 9
  hold it (88.9%, marked not met; at eight seeds 94.0%, Phase 3). The class averages of the three (teacher, student, reference) are kept in
  `keep_cavgs/` beside them for inspection.

**The fast gate on both trees** on the final source, last before the status: Release 15 of 15, Debug 15 of 15
(`R/scratch/logs/p4_final_gate_build.log`, `p4_final_gate_build_debug.log`).

**Files changed in this phase:** this plan (the status line, the file-table row, decisions 3 and 7 (c), section 4,
Phase 3's entry condensed, this entry), the report rewritten as the adoption, `doc/algorithms/class_student.md`,
`doc/implementation_notes/README.md`, and the six documents of the amplitude removal (one deleted).

**Exit.** Decision 7 (c) applied (not met as written, waived for streptavidin by ruling R4): the raw-pick route
adopted; the gates of section 5 pass on the final tree; the report written.
