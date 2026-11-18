# The center student: a real-space network that centres raw picks before the class student

Status, 2026-11-18: planned, revised the same day after the review in section 9 (the response is in section 10).
Hans's rulings R1 to R5 are in (section 3). The margins of decision 14 are new and await Hans's acceptance. The data
sets are named at launch (R6). The plan is not launchable before both. The pose student is left out for now (Hans,
2026-11-18; section 4). Base: `8a690cc19` (clean master; the review read this commit). Supersedes: none. Depends on:
solve2D's class-average centring (validated: yes, the production default); `center_ptcls` (validated: yes,
`doc/algorithms/center_ptcls.md`); the class student on raw picks with the alignment after prediction (validated:
yes, `completed/class_student_raw_picks.md`, adopted 2026-11-16); the students' shared strategy helpers (validated:
yes). Nothing here depends on an unvalidated component.

## In one paragraph

Every measurement on the polar students says the same thing: at good shifts they assign classes and poses well, at
the picking error they do not. The HHCNN is equivariant to rotation about the box centre, so a shift (a phase ramp
spread over every ring) is exactly what it cannot absorb. A real-space convolutional network with stride-1, circularly
padded convolutions is equivariant to whole-pixel circular translations by construction. The center student is such
a network, trained per data set on the shifts of a solve2D teacher, which are already centred on each class average's
support (section 3, decision 1). It reads the particle's box as extracted. It outputs a heatmap over a disk of
candidate centres, plus one "outside" bin. From that it commits an absolute shift, replacing the stored one, or falls
back to the stored shift when the centre is outside the disk or at its edge. By default the class student then
searches its usual disk, under the pick prior, about the new start. A narrower search for confident particles is a
separate option: it is adopted only if it earns it, and only through an explicit record of which particles the
network centred (decision 9). The questions are pre-registered and are answered first in the PyTorch prototype,
through the existing class-student validation, before any Fortran is written:
- Does the network centre unseen particles as well as `center_ptcls` within a declared margin?
- Does the class student from the network's centres, over the usual disk, beat raw picks over the disk?
- Does the narrower search add anything on top?

## 1. Why

- **The students fail on raw picks and succeed on centred ones.** On raw picks the class student's most likely class
  agrees with the full search for 23 to 35% of particles, while its candidate list holds that class for 90 to 96%.
  Dropping centring alone cost 5.6, 15.7 and 19.3 points of mode agreement on streptavidin, bgal and trpv1
  (`completed/class_student_raw_picks_report.md`, Phase 0). The pose student's cone holds its answer for 96% of
  particles at 6.4 Å of shift and for 28% at 12.8 Å (`doc/algorithms/center_ptcls.md`).
- **But the old centred route lost downstream.** The class student on raw picks, aligned over the disk against its
  candidate list, beat the route through `center_ptcls` on every set: by 2.6 to 17.4 points of committed-class
  agreement and 0.75 to 4.80 Å at FRC 0.5. A centring step is therefore judged by what it does to the class student,
  not by how many particles it moves within 5 Å (verdict 13 (c)).
- **The picking error is much larger than the students' tolerance:** a median of 8 to 10 Å, a 90th percentile of 16
  to 22 Å and a 99th of about 40 Å on the real sets.
- **What exists does not close the gap cheaply.** `center_ptcls` needs references and a pass over the data. It puts
  0.45 of trpv1's picks within 3 Å and 0.67 within 5 Å, but only 0.22 to 0.25 and 0.45 on bgal. A rotational average
  centres 0.18 within 3 Å. A power-spectrum student is shift-invariant but learns almost nothing. The amplitude
  channel was rejected (2026-11-16). `polar_student_future_directions.md` section 2 already names the remaining
  candidate: "a translation-equivariant real-space network supervised by converged shifts".
- **Two reference points are measured.** A single pass of bgal's picks against the full run's 50 class averages puts
  0.57 within 3 Å, 0.66 within 5 Å and 0.70 within 6.4 Å. That is a strong measured baseline, not a proven ceiling.
  Two converged solve2D runs agree within 3 Å for 69% (trpv1) and 76% (bgal) of the same particles. That is the
  labels' reproducibility, not an accuracy floor.

The prize: once a solve2D teacher and the network's training have been paid for, prediction needs no references and
no further 2D. It costs a few milliseconds per particle-thread, could later serve the stream's chunks, and gives each
particle a width that may, once calibrated, narrow the search. Its accuracy gain over `center_ptcls` may be modest.

## 2. Contract

- **New programs only, until Phase 4.** `learn_center_student` and `predict_centers` (private) and
  `validate_center_student` (`simple_exec`) are new. Until Phase 4 no existing program, default or record changes.
- **Off state.** Without the new programs, refine2D, refine3D, solve2D, solve3D, `center_ptcls`, the class student
  and the pose student are byte-identical to the base at one thread with a fixed seed (`simulated_workflow` 1jxy and
  6vxx, `test=identity`). The pose student never changes in this plan.
- **The raw-pick route is untouched.** `predict_classes` on a project the center student has not centred behaves
  exactly as at the base, byte for byte. That includes fresh imports, whose records hold `shsig` = 0 (decision 9). The
  default centred route needs no change to `predict_classes` at all.
- **The record.** `predict_centers` writes the shift, `shsig` and `shincr` to `ptcl2D`, and to `ptcl3D` when it has
  rows, and nothing else. The shift is replaced by the network's absolute estimate (decision 5). `shincr` is the
  length of the actual change of the stored shift. `shsig` is the sentinel −1 in the default route, so the next
  search takes the disk under the pick prior about the new start. The project has one owner: parts write part files,
  and the master merges in memory and writes the segments once.
- **Same answer at any thread count.** Training is bit-identical at a fixed thread count, as the HHCNN's is: chunks
  whose size follows from the batch alone, and gradients reduced in chunk order. Prediction is the same at any thread
  count and, with `nparts`, equal byte for byte to shared memory.
- **One copy of everything shared.** The part ranges, the distributed prediction master and the file hash stay in
  `simple_strategy_student_common`. The HHCNN, `polar_student`, `polar_feats` and `src/main/pftc` are not edited.
- **No new command-line parameter** unless ruled. The center student's settings are named constants in
  `src/defs/simple_defs_hhcnn.f90`. Existing keys keep their meaning: `nsample`, `mskdiam`, `nthr`, `nparts`,
  `projfile`, `fname`, `infile`.
- **No new CTest entry.** The center student's workflow checks run inside the existing `class_student` high-level
  entry. Its fixture has planted shifts. `SIMPLE_CTEST_BUDGET` stays unchanged
  (`doc/policies/test_environment_policy.md`: a new entry needs the owner's agreement).
- **Layering.** The commanders are thin `execute(cline)` wrappers. The learn and predict steps are strategy routines,
  `run_<x>(cline)` and `<x>_exec(params, build, ...)`, as the class student's are. The network and the student are
  encapsulated types with `new` and `kill`. Large components are `allocatable` (`doc/policies/compile_time_policy.md`).

## 3. Decisions and rulings

### What the code check found (at `cfbad5a26`, re-read at `8a690cc19`)

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
- **The shift's sign.** A particle image is centred by shifting it by minus its stored shift: the class averager
  calls `shift2Dserial(-shift)` (`simple_classaverager_restore.f90`). The stored shift is therefore the centre's
  displacement from the box's origin, pixel N/2+1 on SIMPLE's even grid, in native pixels.
- **A fresh particle record holds `shsig` = 0, not the sentinel.** Particle parameters initialise to zero
  (`simple_ori.f90`, `pparms`). `shift_needs_disk(0, trs)` is false, so an untouched record would be searched over
  the window. `predict_classes` therefore resets every particle's state to the sentinel before its alignment
  (`simple_strategy_class_student.f90`). `shiftmap_srch%prep` gives the pick prior's weight only to a particle with
  `shsig` < 0, so a finite width also switches the prior off, even when the disk is searched.
- **`record_shift_state` is not a unit conversion.** Its `next_shift_width` inflates the width when the move reached
  the searched radius: a posterior truncated at a bound says nothing of the mass beyond it.
- **The held-out splitter** (`student_holdout_split`) holds out stacks, and falls back to every tenth selected
  particle when there are fewer than ten stacks. It does not guarantee that whole micrographs are held out.
- **The side-artifact pattern exists.** `simple_class_candidates` is a file registered in the project's `out`
  segment and tied to the particle rows by `ptcl_layout_digest`. It is refused when the rows or the project differ,
  and the consumer then runs as without it.
- **`validate_class_student` already has the read at the reference's shifts** (its `predict_dense` step: the network
  at the full search's shifts, no alignment). That is the σ = 0 point of the tolerance curve in Phase 0.
- **The phase flip exists** (`image%apply_ctf(tfun, 'flip', ctfparms)`, used by the class averager).
  `lpstages_setlims` is not the crop policy here: it targets max(2 Å, lp/3) but imposes an 88-pixel floor and FFT
  rounding. Decision 2 sets the center student's own crop.
- **No real-space network exists in SIMPLE.** The HHCNN's machinery is the pattern to copy, not a framework to
  generalise: im2col plus sgemm (team-safe BLAS, explicit loops otherwise), the chunked deterministic reduction, Adam
  and the weights file with a magic number and a version.

### Decisions

1. (2026-11-18, Hans and Claude; R3) **The label is the teacher's final shift, as stored**, from solve2D on `nsample`.
   No correction for each class's offset is needed: solve2D's seg centring already puts the members in the
   class-free support-centroid frame. Training rows are the teacher's particles in classes of more than
   `MINCLSPOPLIM` members (the others were never centred). The label is the centre's displacement from the origin
   pixel N/2+1 of the box as extracted, in Å (native pixels times `smpd`), with the stored shift's sign. Phase 0
   checks the frame on the final class averages (decision 10 (d)).
2. **The input frame.** The particle's box as extracted, with no stored shift applied. It is phase-flipped and
   low-passed at `lp_in` = 12 Å (the pass band of `center_ptcls`, measured best there). It is Fourier-cropped to
   `box_in` = the even number nearest box·smpd/4 Å, so that the realised sampling is `smpd_in` = box·smpd/`box_in`,
   recorded per set. It is normalised to zero mean and unit standard deviation over the background outside
   `mskdiam`. Phase 1 compares a 3 Å target and freezes one. The preprocessing is written once, in Fortran:
   - as Phase 1's exporter, a variant patch that writes the network's input tensors;
   - moved into the tree unchanged in Phase 3.

   So the prototype and the port read identical inputs. The label is converted to the heatmap's pixels with the
   realised `smpd_in`, never with the target 4 Å.
3. **The output: a heatmap over a support disk, plus an outside bin.** One logit per pixel inside the support disk
   of radius `R_s` = min(60 Å, the box's half width − 2 px), about the origin pixel. One more logit means "the centre
   lies outside the support". The softmax runs over all of them. 60 Å is 1.5 times the 99th percentile of the
   picking error (guessed). Phase 0 measures the share of labels beyond `R_s` per set.
   - A training label within `R_s` − 2σ_target gets the Gaussian target of decision 4. A label beyond it gets all its
     target mass on the outside bin.
   - Nothing is renormalised onto the disk, so a centre beyond the support cannot come out as a narrow peak at the
     edge.
   - Regressing (dx, dy) directly was considered and not chosen: under noise it shrinks the estimate toward the pick
     and has no width.
4. **The target** is a normalised Gaussian over the support's pixels about the label, of width
   σ_target = max(σ_label, half a pixel of `smpd_in`). The loss is the cross-entropy of the softmax against it.
   σ_label starts at 1.5 Å, from the reproducibility of two converged solve2D runs. That figure assumes the two runs'
   errors are independent, isotropic, Gaussian, unbiased and of equal variance:
   - P(|d| < 3 Å) = 1 − exp(−9/(2s²)) gives s = 1.96 Å (trpv1, 0.69) and 1.78 Å (bgal, 0.76);
   - a single run's per-axis noise is then s/√2 = 1.39 and 1.26 Å.

   Phase 0 re-measures it from the full difference distribution, not one radius, and freezes it. The target is a
   smoothing of hard labels, not the physical posterior. The width's meaning comes from decision 5's calibration,
   not from the target.
5. **The readout and the commit.**
   - **Estimate.** ŝ is the posterior mean over the 5 × 5 support pixels about the pixel of largest mass, with the
     outside bin left out.
   - **Width about the estimate.** σ_c = sqrt(½ Σ_x p(x) |x − ŝ|²) over the support, after one temperature T scales
     all the logits. T is fitted by negative log-likelihood of the yardstick's centres on the calibration cohort
     (decision 8) and frozen with the model. σ_c is the root-mean-square distance of the posterior from the committed
     estimate, not the spread about the posterior's mean. Two modes 12 Å apart with masses 0.51 and 0.49 give
     σ_c ≈ 5.9 Å.
   - **Validity.** A prediction is invalid when the outside bin's mass is at least 0.5, or when ŝ lies within
     2σ_target of the support's edge (guessed, frozen in Phase 1). An invalid prediction leaves the stored shift
     unchanged and writes the sentinel. It is counted as a fallback and is never dropped from a read.
   - **Commit.** Absolute and replacing: the stored shift becomes ŝ in native pixels (divided by `smpd`), and
     `shincr` is the length of the change. The second-pass variant (Phase 1) prepares the input at the first
     estimate ŝ₁ (shifted by −ŝ₁, as SIMPLE centres an image) and commits ŝ₁ + ŝ₂. The network is the same, since a
     re-prepared box is a translated box (decision 6).
   - **Width in the record.** In the default route the record holds the sentinel, and σ_c goes to the prediction
     table only. A width is written only under decision 9's adaptive route.
6. **The network: a dilated stride-1 trunk with circular padding** (guessed sizes; Phase 1 may change them within 2×
   and freezes them).
   - Seven 3 × 3 convolutions of 16 channels with ReLU, dilations 1, 2, 4, 8, 16, 1, 1, and a 1 × 1 head.
   - Receptive field 1 + 2·(1+2+4+8+16+1+1) = 67 px, about 270 Å at 4 Å.
   - About 1.4 × 10⁴ weights and about 6 × 10⁷ multiply-adds per particle at `box_in` 64.
   - With circular padding the trunk is exactly equivariant to whole-pixel circular translations: a translated input
     gives the translated logit map, at any position.
   - The fixed support and outside bin break equivariance only at the support's edge. That edge is tested
     separately, as the validity rule.
   - Real kernels, and the same deterministic im2col products and chunked reduction as `simple_hhcnn`.
   - A pooled encoder-decoder (stride-2 pooling with upsampling, cheaper but only approximately equivariant) is a
     Phase 1 variant. If it is chosen, its equivariance is measured over every pooling phase and over shifts of a
     quarter and half a pixel, and frozen as a tolerance.
7. **Augmentation.** The eight elements of D4, exact on the pixel grid and unbiased since the in-plane angles are
   uniform. Each is an affine index map about the origin pixel N/2+1 with circular wrap (i → N+2−i mod N for a
   reflection), applied to the image. The label gets the same linear map. A plain array reversal would move the
   origin by one pixel and is not used. There is no translation augmentation in the base design: the real picks'
   displacements supply it. A circular translation over 10 Å with a cosine-tapered box edge (so that the seam is
   smooth) is a Phase 1 variant.
8. **Cohorts, frozen in Phase 0 before any training.** Per set, the micrographs are permuted by a keyed seed and
   split by micrograph identity: T 60%, S 20%, E 20% (guessed).
   - The teacher's `nsample` is drawn from T.
   - Early stopping uses `student_holdout_split` inside the teacher's particles. Its actual grouping (stacks, or the
     every-tenth-particle fallback) is recorded.
   - S is the selection and calibration cohort: Phase 1's variants are chosen on S, and T is fitted on S.
   - E is touched once per frozen design, for the verdicts.
   - Phase 0 checks that the project's stacks map one-to-one to micrographs and records the identity used. If they do
     not, the micrograph name is used.
   - Training: Adam (β 0.9 and 0.999, ε 10⁻⁸, as the HHCNN), batches of 256, a keyed seed. Training stops when the
     held-out loss has not improved for a fixed number of epochs (Phase 1 freezes it).
9. **The consumer.**
   - **Default route (centred + disk).** `predict_centers` then `predict_classes`, unchanged: the reset sends every
     particle to the disk under the pick prior, about its new start. Phase 4 adds this arm to
     `validate_class_student` and changes nothing in `predict_classes`.
   - **Adaptive route, only if verdict 13 (d) keeps it.** Particles whose width sends them to the window
     (3·σ_c ≤ `trs`) keep their width; all others get the sentinel. The record alone cannot say which particles the
     network centred, because fresh records hold 0. So `predict_centers` also publishes a side artifact
     `center_state`, on the `class_candidates` pattern: registered in the `out` segment, tied by `ptcl_layout_digest`,
     listing the particles it centred and their widths. `predict_classes` keeps a stored width only for particles
     listed in a valid artifact, and resets every other particle as today. A missing, stale or mismatched artifact
     means today's behaviour.
   - **Adaptive route's tests:**
     - a direct `predict_classes` on a fresh import (byte-identical to the base);
     - a project centred by `predict_centers`;
     - a legacy record with a zero width and no artifact (reset);
     - a mixed population;
     - an artifact whose digest no longer matches (refused).
10. **Phase 0's measurements, all on the base or a variant patch outside the tree:**
    - (a) per set: box, `smpd`, `mskdiam`, `box_in`, the realised `smpd_in`, and the margin (box·smpd − mskdiam)/2
      against the 99th percentile of the picking error; the share of yardstick centres beyond `R_s`;
    - (b) the yardstick: a full solve2D on every active particle at `SIMPLE_SEED` 1 and 2. Their full difference
      distribution after the common offset (the mean difference vector): its shares within 3 and 5 Å, its
      anisotropy, and σ_label (decision 4);
    - (c) the arms on cohorts S and E against yardstick 1 after the common offset:
      - raw picks;
      - `center_ptcls` (C);
      - the seg operator applied to each particle alone (B0): `calc_cavg_offset` on the phase-flipped particle at
        `lp_in`, with no learning;
      - the strong baseline (O): one pass against yardstick 1's class averages over the disk;
    - (d) the frame: `calc_cavg_offset` on the final class averages of the teacher and of yardstick 1; the share of
      classes under `CENTHRESH`, and the share of particles in classes of 5 members or fewer;
    - (e) the class student's shift tolerance: in `validate_class_student`'s `predict_dense` step, the reference's
      shifts plus isotropic Gaussian noise of σ ∈ {0, 1, 2, 3, 4, 6, 8} Å, keyed by seed. Reads: the mode's agreement
      with the reference, the candidate list holding it, and the aligned class's agreement over the window and over
      the disk;
    - (f) the cohorts of decision 8, the micrograph identity, and the splitter's actual grouping.
11. **The teachers and seeds.** solve2D on `nsample` random active particles of cohort T. Eight seeds per set, paired
    by seed (the teacher and the training both change with the seed; the cohorts do not). The data sets are named at
    launch (R6). Candidates are the raw-pick run's three, with its `nsample`:

    | set | `nsample ncls mskdiam` | import project |
    | --- | --- | --- |
    | streptavidin | `5000 30 85` | `/data/10335-streptav/tst/1_import_particles/tst.simple` |
    | bgal | `2500 50 180` | `/data2/bgal_addon/1_import_particles/bgal_addon.simple` |
    | trpv1 | `2500 50 175` | `/data2/trpv1_c2/1_import_particles/trpv1_c2.simple` |

12. **The reads.** Centring reads on a cohort are against yardstick 1, after the common offset between the teacher's
    frame and the yardstick's (the mean difference over the teacher's own particles):
    - the median residual, the shares within 3, 5 and 6.4 Å, and the 90th percentile, with every particle counted:
      a fallback is scored at its stored (raw) position, never dropped and never scored at a clipped label;
    - the share of fallbacks, and their raw residuals;
    - calibration: the reliability curve, meaning the share of residuals within k·σ_c for k = 1, 2 and 3 against the
      share the heatmaps themselves predict. Reported, never compared with a Gaussian's 0.865;
    - the window escape rate: of the particles that 3·σ_c ≤ `trs` would send to the window, the share whose residual
      exceeds the window's radius (`trs`, in the same unit). Overall, and split by low signal (the lowest decile of the input's in-mask variance), several
      modes (a second local maximum holding at least 0.2 of the mass), anisotropy, and nearness to the support's
      edge;
    - the cost per particle (training and prediction separately, thread-seconds and wall time).

    The downstream reads are those of `validate_class_student`, on cohort E: the committed class's agreement with the
    reference, the class averages' FRC 0.5 crossing, and the route's cost.
13. **Verdicts, written now.** Paired over the eight seeds, on cohort E. LB and UB are the one-sided 95% bounds of the
    paired mean difference (t with 7 degrees of freedom, 1.895). The margins are decision 14's. Each verdict's
    branches are evaluated in the order written, and the first that holds applies.
    - (a) **Early stop after Phase 0:** if LB(B0 − C) > −δ5 within 5 Å on every set, the seg operator alone is not
      worse than `center_ptcls` by more than the margin. Stop with the table: the maintainer decides whether a network
      is still wanted. Otherwise go on.
    - (b) **Centring, after Phase 1:**
      - go if, on every set, LB(N − C) > −δ5 within 5 Å, and on at least one set LB(N − B0) > 0;
      - otherwise negative, and the plan closes in `rejected/`, if UB(N − B0) < δ5 on every set;
      - otherwise inconclusive: stop with the table.
    - (c) **The value for the class student, after Phase 1, before any port.** The prototype's centres, exported,
      enter a variant of `validate_class_student` on the same teacher and the same particles, in three arms:
      - A_raw: raw picks over the disk, today's route;
      - A_cd: the centres over the disk, sentinel, pick prior;
      - A_ca: the centres with the adaptive routing of decision 9.

      The port goes ahead if, on every set, LB(A_cd − A_raw) > −δ_agree on committed-class agreement and
      UB(FRC₀.₅(A_cd) − FRC₀.₅(A_raw)) < δ_FRC (Å, lower is better), and on at least one set LB(A_cd − A_raw) > 0 on
      agreement. Otherwise, if A_ca passes those same conditions against A_raw, stop: the gain would come from the
      search policy, not the centres, and the maintainer decides. Otherwise negative: the plan closes in `rejected/`.
    - (d) **The adaptive route, decided on Phase 1's arms.** It is kept for Phase 4 if, on every set:
      - A_ca is non-inferior to A_cd on agreement and FRC by the same margins;
      - the escape rate's one-sided 95% upper bound is at most ε.

      Otherwise the route is centred + disk and decision 9's artifact is not built. Costs are reported, never gated
      (Hans, 2026-11-15).
    - (e) **The Fortran port, after Phases 2 and 3.**
      - Implementation parity first. On a frozen prototype model and frozen inputs, these must agree with the
        prototype within decision 14's tolerances: preprocessing, the trunk's activations at every layer, logits,
        probabilities, readouts (ŝ, σ_c, validity), loss, gradients, and one Adam update. Any miss is a stop.
      - Then the Fortran-trained models against the prototype's: LB(Fortran − prototype) > −δ5 within 5 Å on every
        set, and the same downstream margins on A_cd. Otherwise stop.
    - (f) **The route in the tree, after Phase 4:** `validate_class_student`'s centred arm (A_cd, or A_ca if (d)
      kept it) against its raw arm, by (c)'s rule. Adopted if it passes; otherwise stop with the table.
14. **Margins and tolerances** (guessed; to be accepted by Hans; frozen before Phase 1 and never moved):
    - δ5 = 0.03 (share within 5 Å);
    - δ_agree = 2 percentage points (committed-class agreement);
    - δ_FRC = 0.3 Å (FRC 0.5 crossing);
    - ε = 0.05 (window escape rate);
    - parity, all in single precision:
      - preprocessing: maximum absolute difference ≤ 10⁻⁵ × the input's root-mean-square;
      - activations and logits: ≤ 10⁻⁴ × the layer's maximum absolute value;
      - probabilities: ≤ 10⁻⁵ absolute;
      - ŝ and σ_c: ≤ 0.01 Å; validity identical;
      - loss: ≤ 10⁻⁵ relative;
      - each gradient tensor: ≤ 10⁻⁴ of its norm;
      - weights after one Adam step: ≤ 10⁻⁵ relative.

### Rulings (2026-11-18, Hans)

- **R1, the prototype: PyTorch first.** Phase 1 is an off-tree PyTorch prototype on exported tensors, outside the
  repository (`~/src/agent_runs/<run>/scratch`). It settles decisions 2, 3, 5, 6 and 8 and gives verdicts 13 (b),
  (c) and (d) before any Fortran is written.
- **R2, the scope: one model per data set for now.** A general model across data sets is a later plan of its own.
- **R3, the teacher: solve2D on `nsample`** for the class student (decision 1). For the pose student on raw picks
  Hans named the shifts solve2D and solve3D accumulated; that route is left out of this plan (section 4).
- **R4, the integration: the class student's route only** (decision 9, verdicts 13 (c), (d) and (f)). The pose
  student, `center_ptcls`, solve2D's and solve3D's own centring, import's defaults, the stream and the picker stay out
  (section 4).
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
- the artifact's name `center_state`;
- the times file kept separate from the prediction table (the class student's reason: the table must be equal
  between parts and shared memory);
- the readout's 5 × 5 window;
- the cohort shares 60/20/20.

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
- Import's and injection's centring defaults (setting the sentinel at import would widen R4's scope), solve2D's and
  solve3D's own centring, the stream, the picker, re-extraction (Joe's item) and `center_ptcls` itself. Replacing any
  of them is a later plan, once verdict 13 (f) is in.
- A general model across data sets, simulated pre-training, and continuous rotation augmentation (R2's later plan).
- Using the heatmap's entropy as a junk score for `planned/selection_student.md`: reported only, not wired.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds Release and Debug, fast gate, the base's identity twice, decision 10's tables and the cohorts frozen, decision 14 accepted, verdict 13 (a) |
| 1 | experiment, off-tree | the exporter's planted-shift and D4 checks; the prototype's reads on S (selection) and E (verdicts) for every set and seed; decisions 2, 3, 5, 6 and 8 frozen; verdicts 13 (b), (c) and (d) |
| 2 | identical (new code, no caller) | the network's tester with plants; layer-wise parity against the frozen prototype (13 (e), first part); fast gate; `unit_numerics`; off-state identity |
| 3 | identical for every existing program; new programs | the student's tester with plants; preprocessing and readout parity; the center student's checks inside the `class_student` entry; parts equal to shared memory byte for byte; training bit-identical at a fixed thread count; verdict 13 (e); off-state identity |
| 4 | result-changing for `validate_class_student` (and for `predict_classes` only if 13 (d) kept the adaptive route) | the raw-pick route byte-identical to the base on a fresh import; decision 9's tests if the artifact is built; `class_student` and `pose_student` workflow tests; the pose student's model byte-identical; refine2D and refine3D identity; verdict 13 (f) |
| 5 | adoption or record | Phase 4's gates on the adopted tree; the report |

## 6. Acceptance

- **Floors and where they come from.**
  - C's, B0's, O's and raw picks' reads on S and E: measured in Phase 0 on the base, frozen in Progress before
    Phase 1.
  - σ_label: measured in Phase 0 (decision 4's derivation from the docs is the starting value only).
  - The class student's tolerance curve: measured in Phase 0, reported.
  - The margins and parity tolerances: decision 14, guessed, accepted and frozen before Phase 1.
  - The prototype's frozen model, inputs and per-layer outputs: produced in Phase 1, frozen before Phase 2.
  - A_raw's reads: measured in Phase 1 beside A_cd and A_ca on the same teacher, and again in Phase 4 beside the
    in-tree arm.
  - The cost estimate (6 × 10⁷ multiply-adds, a few milliseconds per particle-thread, against about 27 ms for
    `center_ptcls`: 40 s × 24 threads over 35,645 particles): guessed, reported, never gated.
- **Identity.**
  - `simulated_workflow` 1jxy and 6vxx at one thread, `SIMPLE_SEED` 7, `test=identity`, base against base first.
    Exclusions only from that comparison.
  - `predict_classes` on a fresh import: byte-identical to the base through the whole plan.
  - The class student's model is byte-identical to the base on bgal at seed 1 (header and weights) through Phase 3;
    the pose student's through the whole plan.
- **Paired bounds** come from the eight seeds (measured). Fewer than eight on a set makes its result "a lead, not a
  finding".
- **Reported, not gated:**
  - the class averages of A_raw and the adopted centred arm at seed 1 for every set
    (`<run>/scratch/keep/validate/<arm>/<set>_s1/`);
  - montages of 64 cohort-E particles with the committed centre, the yardstick's centre and the heatmap, per set,
    including fallbacks and multi-mode cases (the maintainer's eye);
  - the heatmaps' entropy against the class-average rejection's verdict, where the set has one.

## 7. Common rules and stop conditions

- No new parameter unless ruled. Experimental variants are patches on a `git archive` export, not in the tree:
  Phase 0's B0 and tolerance curve, Phase 1's exporter, and its downstream arms. No floor, margin or tolerance
  relaxed without a ruling. A file outside the table gets a row before the edit. The thread policy and the disk policy
  of `how2run_autonomous_agents.md`.
- **Stops:**
  - verdict 13 (a), or an inconclusive or stop branch of (b), (c), (e) or (f), ends at the maintainer;
  - an edit needed in the HHCNN, `polar_student`, `polar_feats`, `src/main/pftc` or the pose student;
  - a new CTest entry needed (an owner decision on the process budget);
  - `predict_classes` on a fresh import no longer byte-identical to the base;
  - training that cannot be made bit-identical at a fixed thread count;
  - a margin (decision 10 (a)) below the picking error's 90th percentile on a set: its box cannot hold the centre and
    needs re-extraction (out of scope);
  - the frame check (decision 10 (d)) finding more than 10% of classes off centre by more than `CENTHRESH`:
    decision 1 would be wrong;
  - cohort E read before a design is frozen;
  - the launch names no data sets (R6), decision 14 is not accepted, or a data set is missing or its checksum changed.
- **Gates no ruling may waive:**
  - the off-state identity;
  - `predict_classes` on a fresh import byte-identical to the base;
  - the class student's identity through Phase 3, and the pose student's throughout;
  - training bit-identical at a fixed thread count;
  - `predict_centers` with parts equal to shared memory byte for byte;
  - parity against the frozen prototype (13 (e), first part).

## 8. Phases

- **Phase 0, the baseline.** Pin the base. Builds Release and Debug, fast gate, the base's identity twice. The data
  sets named at launch (R6), frozen with checksums. Then:
  - the cohorts of decision 8 (micrograph identity checked);
  - decision 10 (a) to (f) on every set at seeds 1 to 8 (yardstick at seeds 1 and 2 only). B0 and the tolerance curve
    run from variant patches;
  - decision 14 accepted by Hans;
  - verdict 13 (a).

  Exit: decision 10's tables, the cohorts and σ_label frozen in Progress; verdict 13 (a) recorded.
- **Phase 1, the prototype** (R1).
  - **The exporter.** A variant patch implementing decision 2's preprocessing in Fortran. It writes each set's input
    tensors, labels (converted with the realised `smpd_in`), yardstick centres and cohort tags. Checked before any
    training:
    - a particle with a planted shift comes out with its label at the right heatmap pixel, sign included;
    - the D4 index maps move the image and the label together (the origin pixel N/2+1 fixed);
    - a non-integral `box`/`box_in` ratio converts exactly.
  - **The training.** PyTorch trains decision 6's network with decisions 3, 4, 5, 7 and 8 at the eight seeds.
    Variants are judged on S, one change at a time:
    - a 3 Å target sampling against 4 Å;
    - the channels halved or doubled;
    - the pooled encoder-decoder (with its equivariance measured);
    - the tapered translation augmentation;
    - the second pass;
    - `nsample` doubled on the largest set.
  - **Freeze.** The chosen design, T and the validity thresholds are frozen. Then the frozen design is read once on E
    for verdict 13 (b).
  - **Downstream.** The frozen design's centres for E enter a variant of `validate_class_student` for verdicts 13 (c)
    and (d).
  - **For the port.** One trained model per set is frozen with its inputs and every layer's outputs for 13 (e).

  Exit: the frozen design and its artifacts in Progress; verdicts 13 (b), (c) and (d).
- **Phase 2, the network.** `simple_rscnn`: real dilated 3 × 3 convolutions with circular padding, ReLU, the 1 × 1
  head, the softmax over the support and the outside bin, cross-entropy against a target map, temperature scaling,
  Adam, the deterministic chunked reduction, and the weights file. If Phase 1 chose the pooled variant, its layers
  instead. Its tester checks:
  - finite-difference gradients;
  - exact equivariance of the trunk to circular whole-pixel translations, separately from the support and outside bin;
  - D4 equivariance of the forward pass with D4-symmetrised weights;
  - bit identity of a training step at 1 and 4 threads;
  - the weights file read back byte for byte;
  - layer-wise parity with the frozen prototype model (decision 14).

  Plants: a wrong gradient sign, a transposed kernel, a padding off by one, and a thread-order reduction. Each must
  be caught. Exit: the tester and its plants, fast gate, `unit_numerics`, off-state identity.
- **Phase 3, the student and its programs.** `simple_center_student`: decision 2's preprocessing (the exporter's
  code, moved), decision 7's augmentation, decision 4's targets, decision 8's training, decision 5's readout, the
  model file (with T and the validity thresholds). Then:
  - the strategy routines `learn_center_student` and `predict_centers` (with parts through `student_predict_parts`),
    and the commanders, UI rows and exec routes;
  - `validate_center_student`: cohorts, teacher, learning, prediction on S and E, yardstick, decision 12's centring
    reads, report and details files;
  - the center student's checks inside the existing `class_student` high-level entry: the model reads back, every
    unseen particle gets a prediction or a counted fallback, the residuals against the planted truth are below the
    raw picks' with LB > 0, and parts equal shared memory.

  Exit: section 5's Phase 3 gates and verdict 13 (e).
- **Phase 4, the consumer.**
  - `validate_class_student` gains the centred arm next to the raw arm on one teacher: A_cd, which needs no change to
    `predict_classes`, or A_ca if 13 (d) kept it.
  - Only for A_ca: the `center_state` artifact and decision 9's change and tests in `predict_classes`.
  - The raw-pick route on a fresh import, byte-identical.

  Exit: section 5's Phase 4 gates and verdict 13 (f).
- **Phase 5, adoption or record.** If the centred route is adopted: `doc/algorithms/class_student.md` (the route),
  a new `doc/algorithms/center_student.md`, the test policy, and `polar_student_future_directions.md` section 2
  (closed for 2D; the pose student's case stays open). If not: the record, with the plan moved to `rejected/` in the
  commit that removes the code. Exit: the report `center_student_report.md` beside the plan.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/learn/simple_rscnn.f90`, `simple_rscnn_tester.f90` (new) | the real-space network | 2, 3, 5 |
| `src/main/learn/simple_center_student.f90`, `simple_center_student_tester.f90` (new) | preprocessing, augmentation, targets, training, readout, calibration, model file | 3, 4, 5 |
| `src/main/strategies/parallelization/simple_strategy_center_student.f90` (new) | `learn_center_student`, `predict_centers` and their master; the cohorts | 3, 4, 5 |
| `src/main/commanders/simple/simple_commanders_center_student.f90` (new) | thin wrappers; `validate_center_student` | 3, 4, 5 |
| `src/main/strategies/parallelization/simple_strategy_student_common.f90` | shared helpers, if needed | 3 |
| `src/main/project/simple_center_state.f90`, `simple_center_state_tester.f90` (new, only if 13 (d) keeps the adaptive route) | the `center_state` artifact | 4, 5 |
| `src/main/strategies/parallelization/simple_strategy_class_student.f90` (only if 13 (d) keeps the adaptive route) | decision 9 | 4, 5 |
| `src/main/commanders/simple/simple_commanders_class_student.f90` | the centred arm beside the raw arm, their reads | 4, 5 |
| `src/defs/simple_defs_hhcnn.f90` | `CENTER_STUDENT_*` constants | 2, 3, 4, 5 |
| `src/main/ui/simple_private_prgs.f90`, `src/main/ui/simple/simple_ui_refine2D.f90`, `production/simple_private_exec_driver.f90`, `src/main/exec/simple_exec_refine2D.f90` | programs, UI and routes | 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` | the unit suites | 2, 3, 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_class_student.f90` | the center student's checks inside the `class_student` entry; the route's reads | 3, 4 |
| `doc/policies/test_environment_policy.md` | suites, the `class_student` entry's checks | 2, 3, 4, 5 |
| `doc/algorithms/class_student.md`, `doc/algorithms/center_student.md` (new) | the route | 5 |
| `doc/implementation_notes/README.md`, `doc/implementation_notes/planned/polar_student_future_directions.md` | record | 0 to 5 |
| `doc/implementation_notes/planned/center_student.md`, `center_student_report.md` (new) | this plan, the report | 0 to 5 |

## 9. Review, 2026-11-18

Reviewed against clean master `8a690cc19`, including the current shift-state producer and consumer, particle
preparation, class centring, held-out splitter and the previous raw-pick experiment. The findings below are review
recommendations, not changes to R1 to R6. No implementation or experiment was launched. The per-data-set,
prototype-first scope is sound; the proposal needs the following contract corrections before implementation.

1. **P1: removing the reset does not preserve the raw-pick route (decision 9).** Fresh particle records initialise
   their static fields to zero (`simple_ori.f90`, `pparms`); import does not generally set `shsig=-1`. The current
   `predict_classes_exec` explicitly sets that sentinel, whereas `shift_needs_disk(0, trs)` is false. Removing the
   reset therefore sends an ordinary imported project to the window without any center-student prediction. Keeping
   the reset in the validation's raw arm hides this regression. Preserve raw-pick behaviour unless valid centring
   state is established, with an explicit contract for old/imported records and prediction provenance. Test a direct
   `predict_classes` call on a fresh import, a student-centred project, a legacy zero-width record and a mixed
   population. Import-wide sentinel changes would widen R4's scope and must not be assumed by this plan.

2. **P1: absolute labels and additive commits use different frames (decisions 1, 2 and 5).** Reading the raw image
   and targeting the final stored teacher shift predicts an absolute displacement from the extraction origin.
   Adding that prediction to an existing shift counts the existing displacement twice. Choose one explicit
   contract: predict absolute shifts from raw input and replace the stored shift, or prepare input at the existing
   shift and train the residual `teacher_shift - input_shift`, then add it back. In both cases `shincr` is the
   actual change of the stored shift. Define the heatmap origin, axis order, sign and conversion through the actual
   Fourier-crop sampling. SIMPLE applies minus the stored shift when preparing an image (`simple_matcher_2Dprep`).
   Include nonzero input shifts, repeated calls, the second-pass variant and non-integral sampling ratios in the
   producer/consumer tests. On an even grid SIMPLE's origin is pixel `N/2+1`; ordinary array reversal moves that
   pixel by one, so D4 image and label transforms need the same affine pixel mapping, not just a vector rotation.

3. **P1: the output disk censors exactly the errors whose uncertainty must stay large (decisions 3 to 5 and 7).**
   The proposal itself reports roughly 40 A errors at the 99th percentile, outside its 30 A output radius;
   augmentation can move further labels outside. A Gaussian renormalised on that disk can become narrow at its
   boundary despite a large true residual. An independent calculation on the proposed 4 A grid, with target
   sigma 2 A and label `(40,0)` A, gives mean `(28,0)` A and `sqrt(trace(Cov)/2)=1.31` A: approximately 12 A error
   with an apparently precise width. Define support for real and augmented labels, an out-of-support policy, and
   a boundary/invalid-prediction fallback. `record_shift_state` is not merely a unit conversion: its
   `next_shift_width` deliberately inflates the width at a search boundary. The new readout needs equivalent
   protection tied to its own support. Do not silently drop tail particles from evaluation or score clipped labels
   as the truth.

4. **P2: the proposed network is not exactly translation-equivariant (decision 6; Phase 2).** Stride-two average
   pooling followed by nearest-neighbour upsampling already fails one-pixel equivariance away from boundaries;
   additive skips do not guarantee its restoration. Four pools retain, at best, a coarser translation subgroup
   before boundary effects. The fixed output disk also prevents unrestricted equivariance of the normalised
   heatmap. Either use an architecture with the required property or describe approximate equivariance and measure
   centre/width stability over every pooling phase and representative subpixel shifts. Separate a trunk test from
   masking/readout tests. Also freeze padding and skip alignment for boxes not divisible by 16. Zero-mean background
   does not remove the seam introduced by a circular translation; test the boundary convention and particle margin.

5. **P2: the width is not yet a calibrated uncertainty about the committed centre (decisions 4, 5, 8 and 12).**
   Cross-entropy against Gaussian-smoothed hard shifts learns a smoothed label distribution, not automatically the
   physical posterior. The label-noise estimate assumes independent, equal-variance, unbiased Gaussian run errors;
   agreement at one radius does not establish those assumptions. Moreover, covariance about the whole posterior's
   mean is not the squared error about the committed local-mode mean. The 0.865 coverage target applies to an
   isotropic Gaussian centred at the estimate, not an arbitrary multimodal heatmap. For masses 0.51 and 0.49 at
   centres 12 A apart, this readout commits the first mode, reports sigma 4.24 A, and contains only 0.51 of the mass
   within twice that sigma. Specify the calibration operation and freeze it on calibration data. Measure actual
   search-window escape rates conditional on being sent to the window, including anisotropic, multimodal, boundary
   and low-signal cases. Distinguish report-only heatmap coverage from evidence sufficient to enable window routing.
   A finite stored width also disables the pick prior even when the disk is selected (`shiftmap_srch%prep`), so
   Phase 4 changes both the starting position and search policy, not just disk versus window.

6. **P2: the verdicts do not establish non-inferiority (decision 13).** Defining worse as a loss beyond two SE makes
   a noisy or underpowered comparison count as "not worse"; for example a four-point loss with SE three points
   passes that rule. Predeclare an acceptable loss and require its one-sided confidence bound to satisfy that
   margin, or label these outcomes inconclusive rather than equivalent. Make the go/negative/stop branches mutually
   exclusive: as written, a network indistinguishable from both C and B0 can satisfy both the go and negative
   clauses of 13(b). This recommendation leaves costs report-only, as ruled.

7. **P2: an aggregate two-SE comparison cannot validate the Fortran port (13(c)).** A wrong transpose, padding,
   target transform or optimiser can stay inside between-seed performance variation. Freeze a prototype model and
   inputs and compare preprocessing, intermediate activations/logits, probabilities, readouts, loss, gradients and
   an Adam update against the same weights with justified floating-point tolerances. Then compare independently
   trained models statistically. Seed spread is an efficacy measure, not an implementation-parity tolerance.

8. **P2: freeze selection and evaluation cohorts before tuning (decisions 8, 11; Phase 1).** The proposal judges
   variants by reads on the unseen particles but does not reserve a final cohort untouched by variant selection.
   Use separate selection/calibration and final evaluation roles, or label the prototype reads exploratory and
   confirm the selected design on untouched data. The existing helper holds out stacks, and falls back to every
   tenth selected particle when fewer than ten stacks occur; it does not guarantee a micrograph holdout. Record
   the actual grouping and fallback, and use micrograph identities if micrograph generalisation is claimed.

Before paying for a Fortran port, the off-tree prototype should also export its predictions back to the existing
class-student validation. Compare raw+disk, centred+disk and centred+width-selected search on the same teacher and
particles. The earlier experiment already found a centred route worse downstream, so matching `center_ptcls`
within 5 A alone does not establish value for the stated consumer. The extra centred+disk arm separates the
centres' effect from the width/search-policy effect. Keep this within R4's class-student scope.

Other corrections before launch: the oracle pass is a strong measured baseline, not a proven information ceiling;
two-run agreement is reproducibility, not an accuracy floor. "Before any 2D" and "without references" apply to
prediction after a solve2D teacher and training have already been paid for. `lpstages_setlims` targets 4 A at 12 A
but also imposes a minimum 88-pixel working box and FFT-size rounding; it does not yield exactly 4 A or a 48-to-64
pixel box in general. Specify the CNN's crop policy and record its realised sampling. The new high-level test also
needs its `simple_test_registry.f90` row, CTest registration and corresponding process-budget update, missing from
the file table; "the CTest list (if needed)" is insufficient for a new high-level entry.

Validation of this review: source tracing and independent standard-library Python calculations for the truncated
Gaussian, two-mode coverage, pooling translation counterexample and even-grid reflection. These are design probes,
not compiled SIMPLE tests or measured model performance. No build, fast gate, training, real-data or Linux/BOX run
was performed. All eight findings and their implementation/experimental checks remain open.

## 10. Response to the review (2026-11-18, Claude)

Every finding is accepted, and the decisions above are rewritten in place (their numbers kept). Where each now
stands:

1. **P1, the reset (decision 9).** The default route no longer touches `predict_classes`: centred + disk keeps the
   reset, so every particle gets the disk under the pick prior about its new start. Narrowing the search is a
   separate, earned option (verdict 13 (d)). It reads the width only through a `center_state` artifact tied to the
   particle rows, never from a bare record. A fresh import stays byte-identical to the base: this is a contract item,
   a stop and a gate no ruling may waive. The five tests named in the review are decision 9's. The import defaults
   stay out of scope (section 4).
2. **P1, the frames (decisions 1, 2 and 5).** Absolute prediction from the box as extracted, replacing the stored
   shift. `shincr` is the actual change. The origin is pixel N/2+1, and the sign is the stored shift's (an image is
   centred by minus its shift). The conversion uses the realised `smpd_in`. The second pass re-prepares the input at
   the first estimate and commits the sum. The D4 maps are affine about the origin with circular wrap. Phase 1's
   exporter checks a planted shift, the D4 maps and a non-integral crop before any training.
3. **P1, the support (decisions 3 and 5).** The support is 60 Å, 1.5 times the 99th percentile, capped by the box.
   An outside bin receives the target mass of labels beyond the support, so nothing is renormalised onto the disk. A
   validity rule (outside mass at least 0.5, or an estimate within 2σ_target of the edge) falls back to the stored
   shift and the sentinel. This is the counterpart of `next_shift_width`'s inflation at a bound. Fallbacks are
   counted and scored at their raw position, never dropped. Phase 0 reports the share of labels beyond the support.
4. **P2, equivariance (decision 6, Phase 2).** The base design is now a dilated stride-1 trunk with circular padding,
   exactly equivariant to circular whole-pixel translations, at about four times the pooled network's multiply-adds.
   The support and outside bin are tested apart from the trunk. Translation augmentation is dropped from the base;
   the tapered variant is Phase 1's. A pooled encoder-decoder is a variant only, with its equivariance measured over
   every pooling phase and over shifts of a quarter and half a pixel.
5. **P2, the width (decisions 4, 5 and 12).** The width is now taken about the committed estimate (the two-mode
   example gives 5.9 Å, not 4.2) and calibrated by one temperature fitted on cohort S. It is read as a reliability
   curve, not against a Gaussian's 0.865. Enabling the window needs the escape rate's upper bound to be at most ε,
   split by low signal, several modes, anisotropy and edge. The default route writes the sentinel, so the pick prior
   stays on. The adaptive route writes a width only for particles it sends to the window.
6. **P2, non-inferiority (decisions 13 and 14).** The margins are declared (δ5, δ_agree, δ_FRC, ε) and judged by
   one-sided 95% bounds over the eight paired seeds. Each verdict's branches are ordered, so exactly one applies;
   anything else is inconclusive and stops.
7. **P2, port parity (13 (e)).** Layer-wise parity against a frozen prototype model and frozen inputs:
   preprocessing, activations, logits, probabilities, readouts, loss, gradients and one Adam step, within decision
   14's tolerances. Only then the statistical comparison. The preprocessing is the same Fortran code in both
   (decision 2).
8. **P2, cohorts (decision 8).** Micrograph-level cohorts T, S and E are frozen in Phase 0. Variants are chosen on
   S, the frozen design is read once on E, and reading E early is a stop. The splitter's actual grouping and fallback
   are recorded, and the micrograph identity is checked.

**The experimental recommendation** is now verdict 13 (c), the gate for the port. The prototype's centres go through
the existing class-student validation as raw + disk, centred + disk and centred + adaptive, on the same teacher and
particles. Matching `center_ptcls` within 5 Å (13 (b)) is necessary but not enough. If only the adaptive arm wins,
the run stops, because the gain would then belong to the search policy, not the centres.

**The other corrections:**
- O is called a strong measured baseline, and two-run agreement is called reproducibility (section 1).
- "Before any 2D" and "without references" are qualified as holding after the teacher and the training.
- The crop policy is the center student's own, with the realised sampling recorded, not `lpstages_setlims`.
- The high-level test is folded into the existing `class_student` entry instead of a new CTest entry: the policy
  wants checks added inside existing suites, and a new entry would need the owner's agreement on the process budget.
  So no `simple_test_registry.f90` row or `SIMPLE_CTEST_BUDGET` change is needed. If a separate entry proves
  necessary, that is a stop.

**What this response did not verify.** The source facts it adds were traced at `8a690cc19`:
- `pparms` initialised to 0;
- the class averager's `shift2Dserial(-shift)`;
- the `class_candidates` and `ptcl_layout_digest` pattern;
- the process-budget rule.

The receptive field, the two-mode width and the operation counts are hand calculations. Nothing was built or run.

## Progress

(one entry per phase, written by the run)
