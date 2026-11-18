# The center student: a real-space network that centres raw picks before the class student

Status, 2026-11-18: planned; Hans's rulings R1 to R5 are in (section 3); the data sets are named at launch (R6), and
the plan is not launchable before that. The pose student is left out for now (Hans, 2026-11-18; section 4). Base:
`cfbad5a26` (the working tree also holds an unrelated edit to `planned/lm_optimizer.md` that this plan does not
touch). Supersedes: none. Depends on: solve2D's class-average centring (validated: yes, the production default);
`center_ptcls` (validated: yes, `doc/algorithms/center_ptcls.md`); the class student on raw picks with the alignment
after prediction (validated: yes, `completed/class_student_raw_picks.md`, adopted 2026-11-16); the students' shared
strategy helpers (validated: yes). Nothing here depends on an unvalidated component.

## In one paragraph

Every measurement on the polar students says the same thing: at good shifts they assign classes and poses well,
at the picking error they do not. The HHCNN is equivariant to rotation about the box centre, so a shift (a phase ramp
spread over every ring) is exactly what it cannot absorb. A real-space convolutional network is equivariant to
translation by construction. The center student is such a network, trained per data set on the shifts of a solve2D
teacher, which are already centred on each class average's support (section 3, decision 1). For each picked
particle it outputs a heatmap: the posterior over where the centre lies in the box. The commit is the posterior mean
of its dominant mode and its width, written to the record `center_ptcls` writes (shift, `shsig`, `shincr`). Dynamic
centring and the class student's alignment then search a window for the confident particles and the disk for the
rest. The question is pre-registered: does the network centre unseen particles at least as well as `center_ptcls`,
without references and at a fraction of its cost? And does the class student, fed its shifts, agree with the full
search at least as often as on raw picks with today's disk alignment?

## 1. Why

- **The students fail on raw picks and succeed on centred ones.** On raw picks the class student's most likely class
  agrees with the full search for 23 to 35% of particles, while its candidate list holds that class for 90 to 96%.
  Dropping centring alone cost 5.6, 15.7 and 19.3 points of mode agreement on streptavidin, bgal and trpv1
  (`completed/class_student_raw_picks_report.md`, Phase 0). The pose student's cone holds its answer for 96% of
  particles at 6.4 Å of shift and for 28% at 12.8 Å (`doc/algorithms/center_ptcls.md`).
- **The picking error is much larger than that tolerance:** a median of 8 to 10 Å, a 90th percentile of 16 to 22 Å
  and a 99th of about 40 Å on the real sets.
- **What exists does not close the gap cheaply.** `center_ptcls` needs references and a pass over the data. It puts
  0.45 of trpv1's picks within 3 Å and 0.67 within 5 Å, but only 0.22 to 0.25 and 0.45 on bgal. A rotational average
  centres 0.18 within 3 Å. A power-spectrum student is shift-invariant but learns almost nothing. The amplitude
  channel was rejected (2026-11-16). `polar_student_future_directions.md` section 2 already names the remaining
  candidate: "a translation-equivariant real-space network supervised by converged shifts".
- **The information ceiling is measured.** A single pass of bgal's picks against the full run's 50 class averages puts
  0.57 within 3 Å, 0.66 within 5 Å and 0.70 within 6.4 Å. No reference-free estimator should be expected to beat it.
  The label's own floor is also measured: two converged solve2D runs agree within 3 Å for 69% (trpv1) and 76% (bgal)
  of the same particles.

The prize: a step that runs before any 2D, costs well under a millisecond per particle, gives a calibrated width for
each particle, and could later serve the stream's chunks and the picker. Its accuracy gain over `center_ptcls` may be
modest, because both are bounded by the same ceiling. Its gains in cost, in needing no references and in the width
are not bounded that way.

## 2. Contract

- **New programs only, until Phase 4.** `learn_center_student` and `predict_centers` (private) and
  `validate_center_student` (`simple_exec`) are new. Until Phase 4 no existing program, default or record changes.
- **Off state.** Without the new programs, refine2D, refine3D, solve2D, solve3D, `center_ptcls`, the class student
  and the pose student are byte-identical to the base at one thread with a fixed seed (`simulated_workflow` 1jxy and
  6vxx, `test=identity`). The class student changes only in Phase 4. The pose student never changes in this plan.
- **The record is `center_ptcls`'s.** `predict_centers` writes only the shift, `shsig` and `shincr` to `ptcl2D`, and to
  `ptcl3D` when it has rows. The class, the angle and every other field are left as they were. The project has one
  owner: parts write part files, and the master merges in memory and writes the two segments once.
- **Same answer at any thread count.** Training is bit-identical at a fixed thread count, as the HHCNN's is: chunks
  whose size follows from the batch alone, and gradients reduced in chunk order. Prediction is the same at any thread
  count and, with `nparts`, equal byte for byte to shared memory.
- **One copy of everything shared.** The held-out split, the part ranges, the distributed prediction master and the
  file hash stay in `simple_strategy_student_common`. The HHCNN, `polar_student`, `polar_feats` and `src/main/pftc`
  are not edited.
- **No new command-line parameter** unless ruled. The center student's settings are named constants in
  `src/defs/simple_defs_hhcnn.f90`. Existing keys keep their meaning: `nsample`, `mskdiam`, `nthr`, `nparts`,
  `projfile`, `fname`, `infile`.
- **Layering.** The commanders are thin `execute(cline)` wrappers. The learn and predict steps are strategy routines,
  `run_<x>(cline)` and `<x>_exec(params, build, ...)`, as the class student's are. The network and the student are
  encapsulated types with `new` and `kill`. Large components are `allocatable` (`doc/policies/compile_time_policy.md`).

## 3. Decisions and rulings

### What the code check found (at `cfbad5a26`)

- **The teacher's shifts are already centred on each class's support.** solve2D sets `center=yes` and
  `center_type=seg` (`simple_commanders_solve2D.f90`; `solve2D_chunks` and the stream pool do the same). The
  `parameters` default `mass` is overridden. From refine2D iteration 3 on, every class with more than `MINCLSPOPLIM`
  (5) members is centred at reference preparation. `calc_2Dref_offset` calls `calc_cavg_offset`
  (`simple_strategy2D_utils.f90`), which:
  - band-passes the average from cenlp up to a high-pass of max(1.5 cenlp, the mask diameter in Å);
  - applies the hard mask and zeroes the negative values;
  - thresholds with Otsu and keeps the largest connected component;
  - takes that component's unweighted centroid (`masscen_cc`).

  An offset of at least `CENTHRESH` (0.5 px) is applied to the reference and added to the class's members
  (`prep2Dref`, `add_shift2class`). The carried even/odd sums get it through `cavger_set_center_offset`. solve2D's
  first stage runs with `center=no`. Chunk runs inside `solve2D_chunks` set `center=no`. Classes of 5 or fewer
  members are never centred.
- **The label is therefore the centroid of the thresholded support of the particle's class average**, not an
  intensity-weighted centre of mass. It does not depend on the class beyond the class's own projection, so it can be
  learned, and it is the frame the class student's teachers were trained in.
- **`predict_classes` resets every particle's centring state** to the sentinel before its alignment
  (`simple_strategy_class_student.f90`: `shsig` = -1, `shincr` = 0). A width written by a centring step is discarded
  today. `shiftmap_srch%prep` would otherwise read it: `shift_needs_disk(shsig, trs)` takes the disk when `shsig` < 0
  or 3·`shsig` > `trs`, the window otherwise.
- **`validate_class_student` already has the read at the reference's shifts** (its `predict_dense` step: the network
  at the full search's shifts, no alignment). That is the σ = 0 point of the tolerance curve in Phase 0.
- **The phase flip exists** (`image%apply_ctf(tfun, 'flip', ctfparms)`, used by the class averager). The library's
  working-grid rule (`smpd_target = max(2 Å, lp/3)`, `lpstages_setlims`) gives 4 Å at a 12 Å band.
- **No real-space network exists in SIMPLE.** The HHCNN's machinery is the pattern to copy, not a framework to
  generalise: im2col plus sgemm (team-safe BLAS, explicit loops otherwise), the chunked deterministic reduction, Adam
  and the weights file with a magic number and a version.

### Decisions

1. (2026-11-18, Hans and Claude; R3) **The label is the teacher's final shift, as stored**, from solve2D on `nsample`.
   No correction for each class's offset is needed: solve2D's seg centring already puts the members in the
   class-free support-centroid frame. Training rows are the teacher's particles in classes of more than
   `MINCLSPOPLIM` members (the others were never centred). Phase 0 checks the frame on the final class averages
   (decision 10 (d)).
2. **The input frame.** The particle's native box is read and phase-flipped. It is low-passed at `lp_in` = 12 Å (the
   pass band of `center_ptcls`, measured best there). It is Fourier-cropped to `smpd_in` = 4 Å, by the library rule at
   that band. It is normalised to zero mean and unit standard deviation over the background outside `mskdiam`. The
   whole cropped box is the input, `box_in` pixels square (even). Phase 1 compares `smpd_in` 3 and 4 Å and freezes one.
3. **The output is a heatmap.** One logit per input pixel. A softmax over the pixels inside the disk of radius
   `R_out` = min(30 Å, the box's half width minus 1 px) is the posterior over the centre's position. The disk is
   three widths of the 10 Å pick prior (`CLASS_DENSE_SHIFT_SIGMA`), the disk the class student's alignment searches.
   Regressing (dx, dy) directly was considered and not chosen: under noise it shrinks the estimate toward the pick and
   has no width.
4. **The target** is a normalised Gaussian over the output pixels about the label, of width σ_target = max(σ_label,
   half an input pixel), with σ_label = 1.5 Å (measured, rounded up). It comes from the reproducibility of two
   converged solve2D runs. If their difference is an isotropic 2D Gaussian with per-axis width s, then
   P(|d| < 3 Å) = 1 − exp(−9/(2s²)). That gives s = 1.96 Å from trpv1's 0.69 and 1.78 Å from bgal's 0.76, so one run's
   per-axis noise is s/√2 = 1.39 and 1.26 Å. Phase 0 re-measures it and freezes σ_label. At 4 Å pixels σ_target is
   2 Å. The loss is the cross-entropy of the softmax against the target.
5. **The readout.** The committed shift is the posterior mean over the 5 × 5 output pixels about the pixel of largest
   mass. It is added to the stored shift in native pixels. The width is sqrt(½ trace Cov) of the whole posterior over
   the disk, converted to the unit `record_shift_state` writes `shsig` in (Phase 3 checks that unit against a
   `center_ptcls` record). The move `shincr` is the committed shift's length. A particle with two modes gets a wide
   width and so the disk downstream: this is intended. There is no draw and no analytic prior: the posterior is
   learned from training displacements distributed as real picks are.
6. **The network** (guessed sizes; Phase 1 may change them within 2× and freezes them). A fully convolutional
   encoder-decoder on `box_in` (about 48 to 64 px at 4 Å):
   - four levels of 8, 16, 32 and 32 channels, each two 3 × 3 convolutions with ReLU, then a 2 × 2 average pool;
   - back up by nearest-neighbour upsampling with additive skips, one 3 × 3 convolution per level;
   - a 1 × 1 head;
   - about 5 × 10⁴ weights and about 1.5 × 10⁷ multiply-adds per particle.

   The receptive field must cover the particle: Phase 1 records it against `mskdiam`/`smpd_in`. Real kernels and the
   same deterministic im2col products and chunked reduction as `simple_hhcnn`.
7. **Augmentation.** The eight elements of D4 (rotations by 90° and mirrors), applied to the image and the label
   together. They are exact on the pixel grid and unbiased, since the in-plane angles are uniform. Also a whole-pixel
   circular translation, uniform over the disk of `a_aug` = 10 Å (guessed), with the label moved by it; the
   background is zero-mean after decision 2, so the wrap adds no edge. Continuous rotations are not used: they need
   interpolation, which changes the noise.
8. **Training.** The students' held-out cohort (`student_holdout_split`: every tenth micrograph) for early stopping
   and the width's calibration. Adam (β 0.9 and 0.999, ε 10⁻⁸, as the HHCNN), batches of 256, a keyed seed. Training
   stops when the held-out loss has not improved for a fixed number of epochs (Phase 1 freezes it).
9. **The consumer, Phase 4 only.** `predict_classes` no longer resets the centring state. It reads the stored `shsig`
   and `shincr` as every 2D search does: the sentinel is set at import and on injection. `validate_class_student`
   keeps setting the sentinel for its raw-pick arm, as it does now. A particle centred by `predict_centers` with
   3·`shsig` ≤ `trs` is then aligned over the window, and the others over the disk.
10. **Phase 0's measurements, all on the base or a variant patch outside the tree:**
    - (a) per set: box, `smpd`, `mskdiam`, and the margin (box·smpd − mskdiam)/2 against the 99th percentile of the
      picking error;
    - (b) the yardstick: a full solve2D on every active particle at `SIMPLE_SEED` 1 and 2; their agreement within 3
      and 5 Å after the common offset (the mean difference vector); σ_label (decision 4);
    - (c) the arms on the unseen particles against yardstick 1 after the common offset:
      - raw picks;
      - `center_ptcls` (C);
      - the seg operator applied to each particle alone (B0): `calc_cavg_offset` on the phase-flipped particle at
        `lp_in`, with no learning;
      - the oracle (O): one pass against yardstick 1's class averages over the disk;
    - (d) the frame: `calc_cavg_offset` on the final class averages of the teacher and of yardstick 1; the share of
      classes under `CENTHRESH`, and the share of particles in classes of 5 members or fewer;
    - (e) the class student's shift tolerance: in `validate_class_student`'s `predict_dense` step, the reference's
      shifts plus isotropic Gaussian noise of σ ∈ {0, 1, 2, 3, 4, 6, 8} Å, keyed by seed. Reads: the mode's agreement
      with the reference, the candidate list holding it, and the aligned class's agreement when the alignment
      searches the window instead of the disk.
11. **The teachers and seeds.** solve2D on `nsample` random active particles. Eight seeds per set, paired by seed.
    The unseen particles are the others. The data sets are named at launch (R6). Candidates are the raw-pick run's
    three, with its `nsample`:

    | set | `nsample ncls mskdiam` | import project |
    | --- | --- | --- |
    | streptavidin | `5000 30 85` | `/data/10335-streptav/tst/1_import_particles/tst.simple` |
    | bgal | `2500 50 180` | `/data2/bgal_addon/1_import_particles/bgal_addon.simple` |
    | trpv1 | `2500 50 175` | `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple` |

12. **The reads** on the unseen particles, against yardstick 1 after the common offset between the teacher's frame and
    the yardstick's (the mean difference over the teacher's own particles):
    - the median residual, the shares within 3, 5 and 6.4 Å, and the 90th percentile;
    - the width's calibration: the share of residuals within 2·`shsig` against 1 − e⁻² = 0.865 for a 2D Gaussian
      (reported);
    - the share whose width sends them to the window (3·`shsig` ≤ `trs`), and their residuals;
    - the cost per particle (training and prediction separately, thread-seconds and wall time).
13. **Verdicts, written now** (paired by seed; worse means the paired mean difference lies in the worse direction
    beyond 2 SE):
    - (a) **Early stop after Phase 0:** if B0 is not worse than C within 5 Å on every set, the operator alone does the
      job. Stop with the table: the maintainer decides whether a network is still wanted.
    - (b) **Learnability, after Phase 1 (the prototype):** go if the network's share within 5 Å is not worse than C's
      on every set. Negative, and the plan closes in `rejected/`, if it does not beat B0 within 5 Å beyond 2 SE on
      any set. Otherwise stop with the table.
    - (c) **The Fortran port, after Phase 3:** its reads equal the prototype's frozen reads within the prototype's own
      seed spread (2 SE) on every set. Otherwise stop: the port differs from what was judged.
    - (d) **The route, after Phase 4:** `validate_class_student` measures both routes on one teacher: today's raw-pick
      route (raw picks, prediction, alignment over the disk, `R_raw`), and the centred route (`predict_centers` from a
      center student trained on the same teacher, then prediction, then alignment over the window or disk by the
      width, `R_cen`). `R_cen` is adopted if, on every set:
      - its committed-class agreement with the reference is not worse than `R_raw`'s;
      - its class averages' FRC 0.5 crossing is not worse;
      - and on at least one set its agreement is better beyond 2 SE.

      If it is not worse on any read but better on none, stop with the table and the costs: the maintainer decides,
      since costs are reported, not gated (Hans, 2026-11-15). Otherwise the losing route is removed.

### Rulings (2026-11-18, Hans)

- **R1, the prototype: PyTorch first.** Phase 1 is an off-tree PyTorch prototype on exported stacks, outside the
  repository (`~/src/agent_runs/<run>/scratch`). It settles decisions 2, 3, 6 and 8 and gives verdict 13 (b) before
  any Fortran is written.
- **R2, the scope: one model per data set for now.** A general model across data sets is a later plan of its own.
- **R3, the teacher: solve2D on `nsample`** for the class student (decision 1). For the pose student on raw picks
  Hans named the shifts solve2D and solve3D accumulated; that route is left out of this plan (section 4).
- **R4, the integration: the class student's route only** (Phase 4, decision 9, verdict 13 (d)). The pose student,
  `center_ptcls`, solve2D's and solve3D's own centring, the stream and the picker stay out (section 4).
- **R5, the names: as proposed.** The center student (`center`, as `center_ptcls` spells it in code; "centring" in
  prose). The modules are `simple_rscnn` (the network) and `simple_center_student`. The programs are
  `learn_center_student`, `predict_centers` and `validate_center_student`, the last in the refine2D UI beside
  `validate_class_student`.
- **R6, the data: decided at launch** ("when I implement this later"). The launch names the sets with their import
  projects and `nsample ncls mskdiam`. Decision 11's table holds the candidates. A launch without them stops in
  Phase 0.

Defaults taken without a question:
- the constants' names, with prefix `CENTER_STUDENT_` in `simple_defs_hhcnn.f90`;
- the model file's magic number and version 1;
- the times file kept separate from the prediction table (the class student's reason: the table must be equal
  between parts and shared memory);
- the readout's 5 × 5 window;
- `a_aug` 10 Å.

## 4. Not in this plan

- The HHCNN, `polar_student`, `polar_feats`, the pose student and `src/main/pftc` (out of bounds; an edit there is a
  stop).
- **The pose student on raw picks** (Hans, 2026-11-18: left out for the time being, to be addressed later). What a
  later plan starts from:
  - Hans's ruling on its teacher: the shifts solve2D and solve3D accumulated (`ptcl3D`).
  - The code check at `cfbad5a26`: solve3D runs with `center=no` (`simple_commanders_solve3D.f90`), so its map is
    never re-centred. Its particles start from the 2D shifts, so the accumulated shifts sit at the de novo map's
    origin. refine3D's own centring (`calcrefvolshift_and_mapshifts2ptcls`: `calc_shiftcen` at cenlp, mapped to the
    particles by `map3dshift22d`; fully for c1, along z for Cn, off for other point groups, several states and the
    fractional update) is the existing operation that would put those labels in the map's centre-of-mass frame. A
    constant 3D offset of the origin projects to a different 2D offset in every view.
  - `planned/pose_student_physical_raw_picks.md` owns the pose student's raw-pick route.
- solve2D's and solve3D's own centring, the stream, the picker, re-extraction (Joe's item) and `center_ptcls` itself.
  Replacing any of them is a later plan, once verdict 13 (d) is in.
- A general model across data sets, simulated pre-training, and continuous rotation augmentation (R2's later plan).
- Using the heatmap's entropy as a junk score for `planned/selection_student.md`: reported only, not wired.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds Release and Debug, fast gate, the base's identity twice, decision 10's tables frozen, verdict 13 (a) |
| 1 | experiment, off-tree | the prototype's reads on every set and seed, decisions 2, 3, 6 and 8 frozen, verdict 13 (b) |
| 2 | identical (new code, no caller) | the network's tester with plants; fast gate; `unit_numerics`; off-state identity |
| 3 | identical for every existing program; new programs | the student's tester with plants; the workflow test `center_student`; parts equal to shared memory byte for byte; training bit-identical at a fixed thread count; verdict 13 (c); off-state identity |
| 4 | result-changing for `predict_classes` and `validate_class_student` | testers of the consumer (a stored width taken; a sentinel takes the disk); `class_student` and `pose_student` workflow tests; the pose student's model byte-identical; refine2D and refine3D identity; verdict 13 (d) |
| 5 | adoption or record | Phase 4's gates on the adopted tree; the report |

## 6. Acceptance

- **Floors and where they come from.**
  - C's, B0's, O's and raw picks' reads: measured in Phase 0 on the base, frozen in Progress before Phase 1.
  - σ_label: measured in Phase 0 (decision 4's derivation from the docs is the starting value only).
  - The class student's tolerance curve: measured in Phase 0, reported.
  - The prototype's reads: measured in Phase 1, frozen before Phase 3 (verdict 13 (c)).
  - `R_raw`'s reads: measured in Phase 4 beside `R_cen`, on the same teacher.
  - The calibration target 0.865: closed form, reported, not gated.
  - The cost estimate (1.5 × 10⁷ multiply-adds, about 1 ms per particle per thread, against about 27 ms for
    `center_ptcls`: 40 s × 24 threads over 35,645 particles): guessed, reported, never gated.
- **Identity.**
  - `simulated_workflow` 1jxy and 6vxx at one thread, `SIMPLE_SEED` 7, `test=identity`, base against base first.
    Exclusions only from that comparison.
  - The class student's model is byte-identical to the base on bgal at seed 1 (header and weights) through Phase 3;
    the pose student's through the whole plan.
- **Paired standard errors** come from the eight seeds (measured). Fewer than eight on a set makes its result "a
  lead, not a finding".
- **Reported, not gated:**
  - the class averages of `R_raw` and `R_cen` at seed 1 for every set
    (`<run>/scratch/keep/validate/<arm>/<set>_s1/`);
  - montages of 64 unseen particles with the committed centre and the yardstick's marked, per set (the maintainer's
    eye);
  - the heatmaps' entropy against the class-average rejection's verdict, where the set has one.

## 7. Common rules and stop conditions

- No new parameter unless ruled. Experimental variants (Phase 0's B0 and tolerance curve, any Phase 1 exporter) are
  patches on a `git archive` export, not in the tree. No floor or tolerance relaxed without a ruling. A file outside
  the table gets a row before the edit. The thread policy and the disk policy of `how2run_autonomous_agents.md`.
- **Stops:**
  - verdict 13 (a), (b), (c) or (d) ends at the maintainer;
  - an edit needed in the HHCNN, `polar_student`, `polar_feats`, `src/main/pftc` or the pose student;
  - training that cannot be made bit-identical at a fixed thread count;
  - a margin (decision 10 (a)) below the picking error's 90th percentile on a set: its box cannot hold the centre and
    needs re-extraction (out of scope);
  - the frame check (decision 10 (d)) finding more than 10% of classes off centre by more than `CENTHRESH`:
    decision 1 would be wrong;
  - the launch names no data sets (R6), or a data set is missing or its checksum changed.
- **Gates no ruling may waive:**
  - the off-state identity;
  - the class student's identity through Phase 3, and the pose student's throughout;
  - training bit-identical at a fixed thread count;
  - `predict_centers` with parts equal to shared memory byte for byte.

## 8. Phases

- **Phase 0, the baseline.** Pin the base. Builds Release and Debug, fast gate, the base's identity twice. The data
  sets named at launch (R6), frozen with checksums. Then decision 10 (a) to (e) on every set at seeds 1 to 8
  (yardstick at seeds 1 and 2 only). B0 and the tolerance curve run from variant patches. Verdict 13 (a). Exit:
  decision 10's tables and σ_label frozen in Progress; verdict 13 (a) recorded.
- **Phase 1, the prototype** (R1). An exporter writes each set's particles in decision 2's frame (SIMPLE's own programs
  where they suffice: `ctfops`, `scale`, `print_project_field` or `export_relion`; otherwise a variant patch), with
  the teacher's labels, the yardstick's shifts and the held-out cohort. PyTorch trains decision 6's network with
  decisions 3, 4, 7 and 8 at the eight seeds, and is judged by decision 12's reads. Measured beside the base design,
  one change at a time:
  - `smpd_in` 3 against 4 Å;
  - the channels halved or doubled;
  - a second pass (predict, re-centre the input on the estimate, predict again);
  - `nsample` 10,000 on the largest set (does more teacher data help?).

  Exit: decisions 2, 3, 6 and 8 frozen with the chosen values, the prototype's reads frozen, verdict 13 (b).
- **Phase 2, the network.** `simple_rscnn`: real 3 × 3 convolutions, ReLU, average pooling, nearest upsampling with
  skips, the 1 × 1 head, the softmax over a pixel disk, cross-entropy against a target map, Adam, the deterministic
  chunked reduction, and the weights file. Its tester checks:
  - finite-difference gradients;
  - translation equivariance (an input shifted by whole pixels away from the border gives a shifted map);
  - D4 equivariance of the forward pass with D4-symmetrised weights;
  - bit identity of a training step at 1 and 4 threads;
  - the weights file read back byte for byte.

  Plants: a wrong gradient sign, a broken skip, a thread-order reduction; each must be caught. Exit: the tester and
  its plants, fast gate, `unit_numerics`, off-state identity.
- **Phase 3, the student and its programs.** `simple_center_student`: decision 2's input frame, decision 7's
  augmentation, decision 4's targets, decision 8's training, decision 5's readout, the model file. The strategy
  routines `learn_center_student` and `predict_centers` (with parts through `student_predict_parts`), and the
  commanders, UI rows and exec routes. Then `validate_center_student`: teacher, learning, prediction on the unseen
  particles, yardstick, decision 12's reads, report and details files. The workflow test `center_student` on the
  class student's fixture (simulated bgal with planted shifts): it gates that the model reads back, that every
  unseen particle is predicted, that the residuals against the planted truth are below the raw picks' beyond 2 SE,
  and that parts equal shared memory. Exit: section 5's Phase 3 gates and verdict 13 (c).
- **Phase 4, the consumer.** Decision 9 in `predict_classes`. `validate_class_student` measures `R_raw` and `R_cen` on
  one teacher. Testers and the workflow tests. Exit: section 5's Phase 4 gates and verdict 13 (d); the losing route
  removed.
- **Phase 5, adoption or record.** If `R_cen` is adopted: `doc/algorithms/class_student.md` (the route),
  a new `doc/algorithms/center_student.md`, the test policy, and `polar_student_future_directions.md` section 2
  (closed for 2D; the pose student's case stays open). If not: the record, with the plan moved to `rejected/` in the
  commit that removes the code. Exit: the report `center_student_report.md` beside the plan.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/learn/simple_rscnn.f90`, `simple_rscnn_tester.f90` (new) | the real-space network | 2, 3, 5 |
| `src/main/learn/simple_center_student.f90`, `simple_center_student_tester.f90` (new) | input frame, augmentation, targets, training, readout, model file | 3, 4, 5 |
| `src/main/strategies/parallelization/simple_strategy_center_student.f90` (new) | `learn_center_student`, `predict_centers` and their master | 3, 4, 5 |
| `src/main/commanders/simple/simple_commanders_center_student.f90` (new) | thin wrappers; `validate_center_student` | 3, 4, 5 |
| `src/main/strategies/parallelization/simple_strategy_student_common.f90` | shared helpers, if needed | 3 |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90`, `simple_class_student_align_tester.f90` (if needed) | decision 9 | 4, 5 |
| `src/main/commanders/simple/simple_commanders_class_student.f90` | the two routes, their reads | 4, 5 |
| `src/defs/simple_defs_hhcnn.f90` | `CENTER_STUDENT_*` constants | 2, 3, 4, 5 |
| `src/main/ui/simple_private_prgs.f90`, `src/main/ui/simple/simple_ui_refine2D.f90`, `production/simple_private_exec_driver.f90`, `src/main/exec/simple_exec_refine2D.f90` | programs, UI and routes | 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` | the unit suites | 2, 3 |
| `src/main/commanders/test/simple_commanders_test_highlevel_center_student.f90` (new), `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90`, the CTest list (if needed) | the workflow test | 3, 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` | the route's reads | 4 |
| `doc/policies/test_environment_policy.md` | suites and workflow test | 2, 3, 4, 5 |
| `doc/algorithms/class_student.md`, `doc/algorithms/center_student.md` (new) | the route | 5 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/polar_student_future_directions.md` | record | 0 to 5 |
| `doc/implementation_notes/planned/center_student.md`, `center_student_report.md` (new) | this plan, the report | 0 to 5 |

## Progress

(one entry per phase, written by the run)
