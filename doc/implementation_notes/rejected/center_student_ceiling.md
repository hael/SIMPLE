# The center student's ceiling: is there anything left for a centring network to win?

Status, 2026-11-18: **stopped by Hans after Phase 0 and the first seed of M1; the center student rejected**
(`rejected/center_student.md`). Run `center_ceiling` on the dell. Base: `a26f57437`. Supersedes: none. Decides the fate of `rejected/center_student.md` (not launched). Depends on: the class
student on raw picks with the alignment over its shortlist (validated: yes, `completed/class_student_raw_picks.md`);
the pose student's centring policy (validated: yes, ruling R7 of `planned/pose_student_physical_raw_picks.md`, adopted
in `a26f57437`); `validate_class_student`'s read at the reference's shifts (validated: yes).

## In one paragraph

`rejected/center_student.md` proposes a real-space network that centres raw picks before the class student. Two
cheaper measurements bound what any such network could win, and this run makes them before anything is built.
**M1, the class student's ceiling:** give the class student the full search's own shifts (perfect centres), and the
same shifts blurred by the error a good network would make, then align over the disk as today. Compare the class
averages with today's route on raw picks. **M2, the pose student's link:** under ruling R7 the pose student predicts
at 2D shifts, which in the stream come from the class student's alignment, not from a full solve2D. Compare the
pose student's unseen map at the class student's shifts with the map at solve2D's shifts. If perfect centres do not
improve the class averages and the class student's shifts serve the pose student as well as solve2D's, no centring
network has a consumer, and the center student is rejected.

## 1. Why

- The class student on raw picks with the alignment over its shortlist (SK) already gives class averages finer at
  FRC 0.5 than the full search's own (seed 1: streptavidin 12.50 against 13.05 Å, bgal 13.65 against 14.14 Å, trpv1
  14.68 against 21.13 Å; `completed/class_student_raw_picks_report.md`). Where the aligned class is the reference's,
  the aligned shift is the reference's to rounding. Its shortlist holds the reference's class for 94 to 99%.
- The route through `center_ptcls` lost to raw picks + SK by 0.75 to 4.80 Å at FRC 0.5.
- The pose student needs centred particles: on streptavidin its cone held the reference for 70% of particles at the
  2D shifts against 19% on raw picks, and the unseen map went from 8.93 to 4.90 Å at FSC 0.143 (pose_raw, R7). In the
  stream those 2D shifts come from the class student's alignment. That link has never been measured: the pose run's
  2D shifts came from a full solve2D over the whole set.
- By ruling R5 of the pose run (and the class student's experience) the maps and class averages decide; agreement
  reads are diagnostics.

## 2. Contract

- **No change to the tree's source.** Every arm is a variant patch on a `git archive` export under
  `scratch/keep/variants/`. The only tree edits are documents (section 8's file table).
- No new command-line key in the tree; the variants may add developer keys to their own copies.
- The class student's and the pose student's code are read, not edited. Costs are reported, never gated.

## 3. Decisions

1. **Data.** The three real sets of the earlier runs:

   | set | class student (M1): import project, `nsample ncls mskdiam` | pose student (M2): selection project, settings |
   | --- | --- | --- |
   | streptavidin | `/data/10335-streptav/tst/1_import_particles/tst.simple`, `5000 30 85` | `/data/10335-streptav/tst/3_selection/tst.simple`, as pose_raw |
   | bgal | `/data2/bgal_addon/1_import_particles/bgal_addon.simple`, `2500 50 180` | `/data2/bgal_addon/3_selection/bgal_addon.simple`, as pose_raw |
   | trpv1 | `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple`, `2500 50 175` | `/data2/trpv1_c2/3_selection/trpv1_c2.simple`, as pose_raw |

   Checksums recorded in Phase 0; a changed or missing file is a stop.
2. **M1's arms**, all inside one `validate_class_student` run per set and seed (one teacher, one model, one
   reference, one unseen set), seeds 1 to 8:
   - A_raw: today's route (raw picks, the shortlist, SK over the disk under the pick prior);
   - A_σ for σ ∈ {0, 1, 2, 3, 4} Å: every unseen particle's stored shift set to the full search's shift plus an
     isotropic Gaussian error of σ per axis (keyed by the seed and the particle, the same draw in every arm), then
     the prediction and SK over the disk about that start, sentinel and pick prior as `predict_classes` sets them.
   A_σ is the class student behind a centring network of error σ; A_0 is the perfect centring network. For scale: two
   converged solve2D runs agree within 3 Å for 69 to 76% of particles (about 1.3 Å per axis each), and the best
   measured centring baseline puts 57% within 3 Å (about 2.3 Å per axis).
   *As wired in Phase 0 (variant V1):* the arms run after A_raw is complete (its teacher, model, reference, prediction
   at the reference's shifts and class averages), so A_raw is untouched; each arm starts again from the teacher's
   result project, sets the stored shifts, runs `predict_classes` unchanged in `predict_sig<σ>/` and makes its class
   averages there. The error is one standard normal pair per particle from SIMPLE's keyed generator (the run's
   `SIMPLE_SEED` and the particle's row), scaled by σ in every arm. Each arm's reads go to one table,
   `validate_class_student_arms.txt`, with A_raw's row computed by the same code (equal to the report's).
3. **M1's reads** per arm: the class averages' FRC 0.5 crossing (Å) and mean FRC to 15 Å, as
   `validate_class_student` reports them; diagnostics: the shortlist's hold on the reference's class, the mode's and
   the committed class's agreement, the committed shift's error against the reference's (median, 90th percentile).
   *As read from Phase 0:* the FRC 0.5 crossing that decides is the report's line "Class averages vs the teacher's,
   FRC 0.5 (A)", the median over the compared classes of the crossing interpolated between shells (`frc05i_A` in
   the arms table). As in the report, an arm compares the classes holding at least 10 unseen particles both in its
   own assignment and in the reference's, so the compared set can differ by a class or two between arms (`ncmp`).
   *Clarified by ruling R4 (the maintainer's session, 2026-11-18, before any verdict): the arms are compared over one class set.* The
   compared set moving between arms shifts the medians by amounts near the 0.3 Å margin. On bgal at seed 1 it
   held 36 classes for A_raw and 35 for every A_σ, and the reference's own median moved from 13.51 to 13.32 Å with
   nothing else changed. So per validation, the paired reads use the classes holding at least 10 unseen particles
   in the reference, in A_raw and in every A_σ alike, and their count and the classes dropped are recorded. The read
   that decides V1 is the interpolated FRC 0.5 crossing, the median over that common set; the mean FRC to 15 Å over
   the same set is reported beside it. The crossing in whole shells moves in steps of a shell (1.2 Å on streptavidin)
   and is reported only. The per-arm columns as the report computes them stay beside. To make this possible, V1
   writes every arm's per-class FRC reads to `validate_class_student_arms_classes.txt`. The validations made before
   the ruling kept no per-class numbers and were rerun.
4. **M2's arms**, on the pose run's frozen models: the physical teacher's 27 trained students (nine teacher-and-seed
   pairs per set, `/data2/agent_runs/pose_raw/scratch/keep/val/P3/`) predicting the same unseen particles with the
   alignment over the cone, as the pose run's arm C3SK (its variant `scratch/keep/variants/C3` there):
   - C_2D: the unseen particles at their solve2D shifts (the pose run's C3SK input);
   - C_CS: the unseen particles at the class student's aligned shifts.
   Both arms run on one build: the C3 variant's sources (on `fd9af0add`, where the models were trained) linked with
   `-static-libgfortran` (the heap-corruption fix of `68638be28`, link only), plus the patch that reads C_CS's shifts.
   C_2D is rerun on that build, not taken from the pose run's tables; its numbers are compared with the recorded
   C3SK's as a check.
   *As wired in Phase 0 (variant V2):* the pose run's C3 source tree is copied as it is: `fd9af0add` with that run's
   uncommitted pose-student changes, on which the models were trained, plus its C3 patch. V2's patch adds the
   developer key `infile2`, a file of rows "particle, x, y" in pixels. Given it, `validate_pose_student` replaces the
   unseen particles' input shifts with those rows, and every unseen particle must have one (arm C_CS). Without it the
   run is C3SK unchanged (arm C_2D). By ruling R3 the pose run's teachers and models are read in place, never
   written; their checksums are verified before and after the runs.
5. **M2's class students.** Per set, three class students (seeds 1 to 3), each as `validate_class_student` builds
   it at the current base: a solve2D teacher on `nsample` particles at their picked positions (the 2D shifts zeroed,
   the sentinel set), `learn_class_student`, then `predict_classes` with SK on the pose student's unseen particles at
   their picked positions. The teacher's draw excludes the pose student's unseen particles (Phase 0 checks that the
   unseen half is the same for every pose pair of a set; if it is not, the draw excludes their union). Pose pair
   (teacher t, seed s) uses class student t. The class student's settings are decision 1's M1 settings
   (`ncls`, `mskdiam`) on the selection project.
   *Revised in Phase 0 under ruling R2 (how the class students are drawn and built):*
   - The unseen half is the same for every pose pair of a set. The pose run splits its training and unseen
     particles greedily by class score, with no seed, so no union is needed. The halves hold bgal 1 256, strep 1 821
     and trpv1 30 504 particles, all of them selected.
   - The base programs cannot keep given particles out of `validate_class_student`'s training draw. A by-hand chain
     (`selection` on the teacher's result project) would leave the training particles' classes and the class
     records' states in the project, and those change the references `predict_classes` prepares. So the class
     students run on a third small variant, VCS: the base plus the developer key `infile`, a list of rows. Given
     it, those rows are the unseen particles and the training draw takes min(`nsample`, rest) of the other active
     particles. Everything else is the base's `validate_class_student`.
   - The input is the import project, which holds the same particle in every row as the selection project (checked:
     stack, index, defocus and astigmatism agree on every row). It is given the selection's states with the base
     program `selection`. The teacher's particles then enter exactly as in M1: at their picked positions, with the
     import project's untouched fields (shift 0, no class).
   - `nsample` is decision 1's, capped at what the draw has left: bgal 1 500 (the whole remaining pool, below
     2 500), strep 5 000 (the whole pool), trpv1 2 500 of 5 000. The class student's committed shifts of the unseen
     rows (its prediction table's x and y) are C_CS's `infile2`.
6. **M2's reads** per arm: the unseen map's FSC 0.143 and 0.5 crossings (Å) against the teacher's map, as the pose
   run read them (beside `fsc0143_reference_A`); diagnostics: the cone's hold on the reference direction, the
   committed direction within 6°, the committed shift's error, and the class student's input shifts against the
   solve2D shifts (median and 90th percentile, Å).
7. **Verdicts, written now.** One-sided 95% bounds of the paired mean difference: M1 over eight seeds (t₇ = 1.895),
   M2 over nine pairs (t₈ = 1.860). The margins are the center student plan's: δ_FRC = 0.3 Å for class averages, and
   the same 0.3 Å for the unseen map (δ_map). Lower Å is better throughout. The first branch that holds applies.
   - **V1, the class student's ceiling,** on g_σ = FRC₀.₅(A_raw) − FRC₀.₅(A_σ) (positive: the centred start is better;
     FRC₀.₅ is the interpolated crossing's median over the validation's common class set, ruling R4, decision 3):
     - *no headroom* if UB(g_2) < δ_FRC on every set;
     - *headroom* if LB(g_2) > δ_FRC on at least one set (record g_σ for every σ: the size of the prize);
     - otherwise *inconclusive*.
     σ = 2 Å decides because no measured centring reaches it; A_0 is reported as the absolute ceiling.
   - **V2, the pose student's link,** on h = FSC₀.₁₄₃(C_CS) − FSC₀.₁₄₃(C_2D) (positive: the class student's shifts are
     worse):
     - *no gap* if UB(h) < δ_map on every set;
     - *gap* if LB(h) > δ_map on at least one set;
     - otherwise *inconclusive*.
   - **The outcome:**
     - V1 no headroom and V2 no gap: **the center student is rejected.** Phase 2 records it (section 8).
     - V1 headroom or V2 gap: stop with the tables; the maintainer decides whether the center student is rescoped to
       that consumer.
     - Any inconclusive verdict: stop with the tables.

## 4. Not in this plan

- Building any part of the center student, or any change to the class student, the pose student, `center_ptcls`,
  solve2D or solve3D.
- The simulated truth sets (their FSC reads saturate at Nyquist; the real sets decide here).
- Costs as a criterion (reported only).

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement setup | base Release build and fast gate; data checksums; the variants built (each with the fast gate); M1's variant run once on one set at seed 1 with A_raw equal to an unpatched `validate_class_student` at the same seed (byte-equal report lines for the raw route; both at one thread, because the 24-thread solve2D teacher does not repeat itself); M2's C_2D on one pair equal to the pose run's recorded C3SK for that pair (within its own rerun noise: the identical number expected) |
| 1 | measurement | M1's 24 validations; verdict V1 |
| 2 | measurement and record | M2's class students and 54 predictions; verdict V2; the outcome; documents only |

## 6. Acceptance

- Floors: none carried from earlier runs; every compared number is measured in this run on one build per
  measurement. The margins are the center student plan's decision 14 (guessed there; adopted here unchanged).
- Reported, not gated: costs; montages are not needed.

## 7. Common rules and stop conditions

- The thread and disk policies of `doc/how2s/how2run_autonomous_agents.md`. The maintainer works on the machine: never
  his checkout or its `build/`; record `uptime`'s load with every timing.
- **Stops:** a data set missing or changed; a Phase 0 equality check failing; an inconclusive verdict or a headroom or
  gap verdict (section 3, 7); an edit needed in the tree's source.

## 8. Phases

- **Phase 0, the setup.** Pin the base, build Release, fast gate. Data checksums. The M1 variant (A_σ arms inside
  `validate_class_student`), the M2 variant (C3 plus the shift input, static runtime), the class-student shift
  producer (decision 5). The equality checks of section 5. Exit: the variants and checks in Progress.
- **Phase 1, M1.** Three sets × eight seeds. Tables per set: every arm's reads and the paired g_σ with bounds.
  Exit: V1 in Progress.
- **Phase 2, M2 and the outcome.** Nine class students, then C_2D and C_CS for the 27 pairs. Tables per set; V2;
  the outcome of decision 7.
  - Rejected: append a section "Rejected, 2026-11-18" to `rejected/center_student.md` with both tables and the
    reasoning in a few sentences, move it to `doc/implementation_notes/rejected/center_student.md` (plain `mv`), and
    update its row in `doc/implementation_notes/README.md`.
  - Otherwise: the tables in this plan's Progress and a status line on `rejected/center_student.md` pointing here.
  - The report `center_student_ceiling_report.md` beside this plan, which stays in `planned/`.
  Exit: the report; `git status --porcelain` listed in it.

### File table

| File | What | Phases |
| --- | --- | --- |
| `doc/implementation_notes/rejected/center_student_ceiling.md`, `center_student_ceiling_report.md` (new) | this plan, its Progress and report | 0, 1, 2 |
| `doc/implementation_notes/rejected/center_student.md` → `doc/implementation_notes/rejected/center_student.md` | the record of the outcome (moved only if rejected) | 2 |
| `doc/implementation_notes/README.md` | the plans' rows | 0, 2 |

## Progress

(one entry per phase, written by the run)

### Phase 0, the setup: done 2026-11-18

Run root `R` = `/data2/agent_runs/center_ceiling`. Logs are under `R/scratch/logs/`, kept outputs under
`R/scratch/keep/`, and each kept item is listed in `R/scratch/keep/summary.tsv`. The run journal is
`R/state/journal_phase_0.md`.

**Base, machine, configuration.**
- Base: `a26f57437d72ed3ad89bfb80f5b26c04eac8d43d` (HEAD, unmoved), exported with `git archive` to
  `R/scratch/keep/src_base`.
- Machine: the dell. Oracle Linux 8.10, 24 cores, 62 GB, no GPU. Compilers: gfortran 15.2.1 (gcc-toolset-15) and
  CMake 3.26.5. The load average was 0.41 at the start; the maintainer ran nothing during the phase.
- Every tree is configured the same way:
  `cmake -S <src> -B <tree> -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DNICE=OFF -DUSE_ARCHOPT=ON -DFORTLS=OFF
  -DPython3_EXECUTABLE=/usr/bin/python3.12`, then `make -j20` (V1, VCS and V2 were built side by side at 20, 12
  and 8 jobs), `make install` and the fast gate. V2 also gets `-DCMAKE_EXE_LINKER_FLAGS=-static-libgfortran`.
- Every script points `SIMPLE_PATH` and `PATH` at the tree it means (`tools/env.sh`). The login environment
  points them at the maintainer's build.

**Builds.** Each passed the fast gate (`scripts/run_fast_gate.sh`, 15 of 15 tests):

| tree | what it is | make wall | log |
| --- | --- | --- | --- |
| `build_base` | the base | 86 s | `build_base.*` |
| `variants/V1/build` | M1's variant | 97 s | `variant_V1.*` |
| `variants/VCS/build` | M2's class-student producer | 140 s | `variant_VCS.*` |
| `variants/V2/build` | M2's pose-student variant | 174 s | `variant_V2.*` |

**The three variants.** Each is a copy of a source tree with its patch beside it, under
`R/scratch/keep/variants/<name>/`. Nothing in the repository's sources changed.
- **V1** (`V1.patch`, 194 lines, one file: `simple_commanders_class_student.f90`).
  - After today's route (A_raw) is complete, `validate_class_student` runs the five centred starts A_σ, σ = 0 to
    4 Å, as decision 2 now records.
  - The arms' reads go to `validate_class_student_arms.txt`, one row per arm. The routine that computes the class
    averages' Fourier ring correlation (FRC) against the teacher takes the arm's classes and stack.
  - The arms' wall time is kept out of the report's run-time table.
- **VCS** (`VCS.patch`, 81 lines). `validate_class_student` with a developer key `infile`, a list of rows. Those
  rows become the unseen particles, and the teacher's draw takes the other active particles (decision 5's
  revision).
- **V2** (`V2.patch`, 60 lines on top of the pose run's C3 tree and its `C3.patch`). `validate_pose_student` with a
  developer key `infile2`, rows "particle, x, y" in pixels, which replace the unseen particles' input shifts
  (decision 4).
  - The C3 tree is `fd9af0add` plus the pose run's uncommitted pose-student changes plus C3's patch, not a bare
    `fd9af0add`. It is copied as it is, since the models were trained on it.
  - Linkage: V2's `simple_exec` holds the Fortran runtime's own I/O routines, as the base's does since `68638be28`.
    C3's did not. `libgfortran.so.5` is still loaded for ARPACK and OpenBLAS, as on the base.

**Data (decision 1).**
- MD5 checksums of the six projects and of the 407 particle stacks they reference are in
  `R/scratch/keep/data_check/data.md5` (`tools/data_check.sh`). They were verified again at the end of the phase:
  no change.
- The three selection projects have the checksums the pose run recorded.
- The three import projects' particle records print byte-equal to the class-student run's dumps (`class_raw`).
- The import and selection projects hold the same particle in every row on all three sets. Stack, index in the
  stack, defocus and astigmatism agree on every row (`tools/rows_match.py`).
- The pose run's nine frozen teachers match their checksums before and after this phase's runs
  (`tools/teacher_check.sh`), and so do its 27 trained models (`data_check/p3_models.md5`). Ruling R3, the pose
  run read only, holds.

**M2's class-student producer (decision 5).**
- The pose student's unseen half is the same for all nine teacher-and-seed pairs of a set. This was read from
  the pose run's prediction tables (`tools/m2_rows.py`, log `m2/rows/m2_rows.log`), and the reason is in the
  code: the pose run splits training from unseen particles greedily by class score, with no random draw.
  - Unseen: bgal 1 256, strep 1 821, trpv1 30 504 particles, all selected.
  - Particles left for the class student's teacher: 1 500, 5 000 and 5 000.
- The driver `tools/cs_m2.sh <set> <t>` copies the import project, gives it the selection's states with the base
  program `selection`, and runs VCS's `validate_class_student` with `SIMPLE_SEED=t` and the unseen rows. It then
  writes the class student's committed shifts as `cs_shifts.txt`.
- Smoke on streptavidin, class student 1: 5 000 training and 1 821 unseen particles, 4 min 35 s on a loaded
  machine, 1 821 shifts written.
- Against solve2D's shifts, the class student's shifts differ by a median of 6.5 Å and a 90th percentile of
  33.5 Å; 28% lie within 3 Å (`tools/cs_vs_2d.py`; a diagnostic of decision 6).

**Equality checks (section 5): both pass.**
- **V1's A_raw against the unpatched base.** Both ran `validate_class_student` on bgal at seed 1 and one thread,
  because the 24-thread solve2D teacher does not repeat itself (`m1/eq_base/bgal_s1`, `m1/eq_V1/bgal_s1`).
  - The report and the details file differ only in the cost and timing lines (and in the cost verdict those lines
    decide).
  - The prediction tables (aligned, and at the reference's shifts), the student's weights, and the student's and
    the reference's class-average stacks are byte-identical.
  - V1's arms table repeats the report's A_raw reads: FRC 0.5 13.19 Å, mean FRC to 15 Å 0.825, aligned class equal
    to the reference's for 54.6%.
  - The drawn start errors have the expected spread: median 1.187σ and 90th percentile 2.114σ, against 1.177σ and
    2.146σ for an isotropic Gaussian.
- **V2's C_2D against the pose run's recorded C3SK** (the alignment over the cone at solve2D's shifts), on the
  streptavidin pair teacher 1, seed 1 (`m2/val/C2D/strep_t1_s1/kept`).
  - The report differs only in timings and costs. Every agreement read and both map rows are identical: the unseen
    map's FSC 0.5 is 8.58 Å and FSC 0.143 4.73 Å, against the reference's 4.16 Å.
  - In the prediction table, 16 of 1 821 rows differ in the fifth decimal of the top probability or the
    confidence. All poses, shifts and directions are identical.

**Smokes (not measurements).**
- V1 at 24 threads, streptavidin seed 1 (`m1/smoke_V1/strep_s1`): 1 min 22 s in all, of which the five arms took
  24 s. The class-average stacks of the teacher, A_raw, A_0 and the reference are kept.
- V2's C_CS on the streptavidin pair teacher 1, seed 1 (`m2/val/CCS/strep_t1_s1/kept`): the 1 821 shifts were
  read, and the unseen map reached FSC 0.143 at 7.22 Å against C_2D's 4.73 Å.
- One pair decides nothing: verdict V2 comes from Phase 2's 27 pairs. Phase 1 runs its own 24 validations.
- Costs on a loaded machine are not read.

**Exit items.**
- The base pinned, built and gated: done (above).
- Data checksums: recorded and re-verified.
- V1, V2 and the producer (VCS with `cs_m2.sh`) built, each through the fast gate. The producer and C_CS were
  smoke-run.
- The unseen half checked: the same for every pair of a set.
- Both equality checks: passed.

**Files changed in the repository:** this plan only. Decisions 2 to 5 now record how Phase 0 wired them; decision 5
is revised under ruling R2; section 5's Phase 0 gate states the one-thread condition; and this entry. The README
row of this plan also changed.

### Phase 1, M1: stopped after seed 1 (2026-11-18)

Hans stopped the run at 18:40 on the first seed's numbers ("move the real-space convolutional network to the
rejected folder and start testing more difficult data"); verdicts V1 and V2 were never formed and M2 was not run.
What was measured, one seed per set (not a verdict: the plan's bounds need eight):

The class averages' interpolated FRC 0.5 crossing (Å, lower is better), the median over the classes compared in
every arm (ruling R4: streptavidin 27 classes, bgal 31); A_σ starts each unseen particle at the full search's shift
plus a Gaussian error of σ Å per axis, A_raw at its picked position; every arm then aligns over the disk:

| set | A_raw | A_0 | A_1 | A_2 | A_3 | A_4 | the reference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| streptavidin | 12.63 | 12.52 | 12.55 | 12.64 | 12.87 | 12.71 | 12.84 |
| bgal | 13.15 | 13.00 | 13.04 | 13.02 | 13.02 | 12.85 | 13.40 |
| trpv1 (before R4: each arm's own class set; its rerun was stopped) | 18.60 | 21.02 | | 19.91 | | | |

- Perfect centres (A_0) gain 0.12 Å on streptavidin and 0.15 Å on bgal, half the 0.3 Å margin, and lose 2.4 Å on
  trpv1. At σ = 2 Å, closer than any measured centring, the gain is −0.01 and 0.13 Å.
- The network's most likely class improves with a centred start (streptavidin 21.5 → 34.0%, bgal 35.6 → 57.4% at
  A_0), but the class committed after the alignment over the shortlist does not (streptavidin 28.5 → 26.2%, bgal
  54.3 → 55.9%): the alignment already recovers what centring would give.
- Every arm's class averages are as fine as or finer than the full search's own on streptavidin and bgal.

Logs, tables and seed 1's class averages: `/data2/agent_runs/center_ceiling/scratch/keep/m1/` (`M1/` after R4,
`M1_preR4/` before) and the journal `state/journal_phase_1.md`.

