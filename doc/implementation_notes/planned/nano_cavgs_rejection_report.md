# Report: a nano context for model_cavgs_rejection

Run of 2026-10-27/28 on the Mac, for the plan `nano_cavgs_rejection.md`. Phase 4 (the 3D criterion) was cancelled
by the maintainer's ruling once the agreement criterion had decided.

## The answer

- What was delivered is the `nano` context of `model_cavgs_rejection`: validity gates, a bank of six lattice
  features, and normalization against the project's floor of lattice-free classes.
- No `nano` preset is promoted. The learned model ranks the classes exactly as one feature, `lattice_peak`, does.
  Its learned threshold costs specificity, so on every held-out dataset it falls below the best single-feature
  baseline.
- The decision to use by hand is a threshold on `lattice_peak` of about 2 floor standard deviations above the
  floor. Pooled over the three labeled datasets, 2.09 floor sd reproduces the maintainer's selection with a
  balanced accuracy of 0.941. That number is in sample.
- Added after the run (2026-10-28): the built-in preset `nano` (`quality_model=nano`) is a threshold on the raw
  `lattice_peak`, not on the floor-normalized one. The run never tried a fixed threshold as a baseline; on the
  three labeled tables it is as good or better: `lattice_peak >= 0.93` gives balanced accuracy 0.984, 0.936 and
  0.939 (plateau 0.90 to 0.95; best per project 0.91, 0.94, 0.845), against 0.967, 0.935 and 0.942 for the
  2.09-floor-sd gate. The ratio is normalized by its ring mean, so the threshold carries across datasets, and it
  needs no floor: on 71-17/traj01 after `refine2D_nano` (chunk classes) the mixture split a continuum (separation
  1.2 floor sd), the guard failed, and a floor-relative preset kept every class; the raw threshold rejects 71 of
  605, those without a lattice by eye. The maintainer then retired the floor altogether (ruling of 2026-10-28):
  the nano bank is `lattice_peak` raw plus median/MAD for the other five, and `thres=` moves the preset's
  threshold for a run. The threshold itself was then set by the maintainer at `LAT_PEAK_MIN` = 0.65, for recall:
  at 0.93 and at 0.90 the rejected pile of 71-17/traj01 still held classes with faint fringes that he keeps, and for
  these time series losing useful frames costs more than carrying some junk into 3D. On the labeled sets 0.65
  keeps 98 to 99% of the selected classes (TPR 0.984, 0.988, 0.981) and lets 27 to 56% of the rejected ones
  through (TNR 0.733, 0.444, 0.670); on traj01 it rejects 7 of 605 valid classes beside the 44 empty ones. The
  learned model is still not promoted.
- The stored class correlation, the old by-hand gate, is clearly worse: balanced accuracy 0.756 to 0.890 with a
  threshold learned on the other datasets, 0.860 to 0.937 with the best threshold on each project itself.

Terms. A class average is the aligned mean of a class's frames. A label is the maintainer's keep/reject state of a
class. The floor is the lower component of a two-component Gaussian mixture fitted to the classes' `lattice_peak`:
the classes without a visible lattice. One "floor sd" is that component's standard deviation. AUC is the chance that
a randomly chosen kept class scores above a randomly chosen rejected one. Balanced accuracy (BA) is the mean of the
fraction of kept classes accepted (TPR) and the fraction of rejected classes rejected (TNR). Classes are counted
once each; the two empty classes of each 71-6 set are left out.

## What was learned

1. The cryo-EM features do not transfer as they are. At 0.358 A per pixel in a 43 to 57 A box, every
   nanoparticle class is hard rejected in the chunk, pool and sieve contexts. The 100-40 A band is empty, and
   the foreground geometry, the band-pass feature and the relational feature are constant.
2. The lattice is the evidence. The per-class AUC of the raw features against the labels:

   | feature | 71-17/traj04 | 71-6/traj01 | 71-6/traj04 |
   | --- | --- | --- | --- |
   | `lattice_peak` | 0.990 | 0.975 | 0.983 |
   | `spot_snr` | 0.987 | 0.972 | 0.974 |
   | `spot_count` | 0.985 | 0.640 | 0.873 |
   | `band_presence` | 0.995 | 0.585 | 0.841 |
   | `spot_arc` (lower is better) | 0.322 | 0.255 | 0.244 |
   | `spot_diam` | 0.214 | 0.372 | 0.302 |
   | stored class correlation | 0.965 | 0.907 | 0.908 |
   | population | 0.227 | 0.437 | 0.322 |
   | relational feature (nano band) | 0.559 | 0.404 | 0.522 |

3. With the cryo-EM resolution constants divided by 8 or 16, the local-variance features become strong (AUC 0.95
   to 0.97). They still add nothing once the lattice family is in: the held-out balanced accuracy changes by less
   than 0.01 without them, against 0.08 to 0.11 without the lattice family. By the plan's rule, only the lattice
   family entered the nano bank.
4. The floor must be guarded, and the guard must fit class averages. On 71-6/traj01 the rejected classes form two
   groups, so the floor is wide and lies only 3.7 floor sd below the lattice classes. The plan's guard of 4 sd
   failed there. The fallback normalization put that dataset on another scale, and every lattice model failed on
   it (balanced accuracy 0.50 at AUC 0.975). The guard is now: the floor owns at least 5% of the classes, lies at
   least 3 floor sd below the upper component, and that gap is at least 0.2 in log10 ratio. The last condition
   catches a narrow split of a single population, seen on a synthetic project where every class carries a
   lattice.
5. The learned model is `lattice_peak` with extra steps. Leave one dataset out, the real learner chose
   probability thresholds of 0.10 to 0.20. Its objective keeps nearly every selected class (recall 0.95 to 0.99)
   and pays in specificity (0.71 to 0.85). At its best threshold on each project, the model reaches what
   `lattice_peak` alone reaches.
6. The maintainer's selection sits at about 2 floor sd, not at the 4 sd of the per-frame floor rule. The 4-sd rule
   keeps only 28% to 75% of the selected classes on the 71-6 sets.

## The held-out tables (every baseline)

Each labeled dataset is held out in turn. The model is learned on the other two by `model_cavgs_rejection
quality_mode=learn`. It scores the held-out project by project-backed `evaluate`. A "learned threshold" for a
single-feature baseline is chosen on the two training datasets, by pooled balanced accuracy.

| method | held out | threshold | AUC | TPR | TNR | balanced accuracy |
| --- | --- | --- | --- | --- | --- | --- |
| nano model (learned, its threshold) | 71-17/traj04 | 0.20 | 0.991 | 0.984 | 0.849 | 0.916 |
| nano model, best threshold on the project itself | 71-17/traj04 | 0.778 | 0.991 | 0.984 | 0.987 | 0.985 |
| stored correlation alone (learned threshold) | 71-17/traj04 | 0.206 | 0.965 | 0.978 | 0.802 | 0.890 |
| `lattice_peak` alone, floor sd (learned threshold) | 71-17/traj04 | 1.997 | 0.985 | 0.984 | 0.935 | 0.960 |
| relational feature alone (learned threshold) | 71-17/traj04 | -1.918 | 0.559 | 1.000 | 0.000 | 0.500 |
| stored correlation, best threshold on the project itself | 71-17/traj04 | 0.2233 | 0.965 | 0.924 | 0.950 | 0.937 |
| lattice floor rule, 4 floor sd | 71-17/traj04 | 0.945 | 0.990 | 0.978 | 0.989 | 0.984 |
| nano model (learned, its threshold) | 71-6/traj01 | 0.15 | 0.973 | 0.953 | 0.852 | 0.903 |
| nano model, best threshold on the project itself | 71-6/traj01 | 0.621 | 0.973 | 0.896 | 0.984 | 0.940 |
| stored correlation alone (learned threshold) | 71-6/traj01 | 0.707 | 0.907 | 0.515 | 0.996 | 0.756 |
| `lattice_peak` alone, floor sd (learned threshold) | 71-6/traj01 | 2.736 | 0.975 | 0.838 | 0.992 | 0.915 |
| relational feature alone (learned threshold) | 71-6/traj01 | 0.088 | 0.405 | 0.405 | 0.428 | 0.416 |
| stored correlation, best threshold on the project itself | 71-6/traj01 | 0.2642 | 0.907 | 0.788 | 0.973 | 0.880 |
| lattice floor rule, 4 floor sd | 71-6/traj01 | 1.339 | 0.975 | 0.280 | 1.000 | 0.640 |
| nano model (learned, its threshold) | 71-6/traj04 | 0.10 | 0.983 | 0.989 | 0.711 | 0.850 |
| nano model, best threshold on the project itself | 71-6/traj04 | 0.363 | 0.983 | 0.943 | 0.957 | 0.950 |
| stored correlation alone (learned threshold) | 71-6/traj04 | 0.287 | 0.908 | 0.830 | 0.871 | 0.850 |
| `lattice_peak` alone, floor sd (learned threshold) | 71-6/traj04 | 2.173 | 0.983 | 0.902 | 0.978 | 0.940 |
| relational feature alone (learned threshold) | 71-6/traj04 | 0.604 | 0.522 | 0.303 | 0.718 | 0.510 |
| stored correlation, best threshold on the project itself | 71-6/traj04 | 0.2297 | 0.908 | 0.822 | 0.897 | 0.860 |
| lattice floor rule, 4 floor sd | 71-6/traj04 | 1.113 | 0.983 | 0.746 | 1.000 | 0.873 |

The promotion test of plan section 8 compares the model with the best baseline on each held-out dataset; the
per-project thresholds are excluded, because they are fitted on the held-out data.

| held out | model AUC | best baseline AUC | model BA | best baseline BA | model not below both |
| --- | --- | --- | --- | --- | --- |
| 71-17/traj04 | 0.991 | 0.990 (floor rule, 4 sd) | 0.916 | 0.984 (floor rule, 4 sd) | no |
| 71-6/traj01 | 0.973 | 0.975 (`lattice_peak`, learned threshold) | 0.903 | 0.915 (`lattice_peak`, learned threshold) | no |
| 71-6/traj04 | 0.983 | 0.983 (`lattice_peak`, learned threshold) | 0.850 | 0.940 (`lattice_peak`, learned threshold) | no |

The gate by hand, `lattice_peak` at 2.09 floor sd, pooled over the three datasets (in sample):

| dataset | TPR | TNR | balanced accuracy |
| --- | --- | --- | --- |
| 71-17/traj04 | 0.984 | 0.950 | 0.967 |
| 71-6/traj01 | 0.882 | 0.988 | 0.935 |
| 71-6/traj04 | 0.909 | 0.978 | 0.944 |

The stored correlation's best thresholds on the three projects are 0.2233, 0.2642 and 0.2297. A threshold near
0.23 to 0.26 is the gate on the stored correlation if no lattice measurement is wanted. It is worse in every case.

## The hold-out trajectories

The seven unlabeled trajectories were run through `analysis2D_nano` (defaults, `element=Pt`). Each project was
copied, and the copy was scored with `apply` (no pruning) by the model learned on all three labeled datasets. The
floor columns come from a fit to the `lattice_peak` the program wrote.

| trajectory | classes | hard rejected | accepted by the model | floor weight | floor separation (floor sd) | gap (log10) | guard | at or above 2.09 floor sd | accepted by both | particle diameter (A) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 71-17/traj01 | 649 | 20 | 383 | 0.45 | 2.56 | 0.52 | fail | 329 | 290 | 27.0 |
| 71-17/traj02 | 649 | 29 | 377 | 0.61 | 3.48 | 0.62 | pass | 253 | 253 | 25.5 |
| 71-17/traj03 | 649 | 18 | 331 | 0.36 | 0.56 | 0.24 | fail | 0 | 0 | 30.0 |
| 71-17/traj05 | 649 | 43 | 311 | 0.68 | 5.20 | 0.64 | pass | 214 | 214 | 25.5 |
| 71-17/traj06 | 649 | 33 | 395 | 0.49 | 5.15 | 0.75 | pass | 336 | 336 | 32.0 |
| 71-6/traj02 | 684 | 31 | 322 | 0.59 | 2.88 | 0.58 | fail | 278 | 276 | 23.0 |
| 71-6/traj03 | 684 | 23 | 414 | 0.48 | 6.66 | 0.46 | pass | 331 | 323 | 20.0 |

- Four trajectories have a floor. There the model accepts every class the 2.09-sd gate accepts, and more.
- 71-17/traj03 has no lattice population: the gap is 0.24 and the separation 0.56 sd. Its classes should be
  judged by eye.
- 71-17/traj01 and 71-6/traj02 have a clear lattice population (gaps 0.52 and 0.58), but their floors are wide
  (2.6 and 2.9 sd). There the lattice features fall back to median/MAD, and the model's decision is not
  calibrated.
- "Particle diameter" is the median `spot_diam` over the classes the model accepted, as `apply` records it.
- The stacks for review (selected, rejected, ranked, hard-gate rejects) and the per-class tables are in
  `~/nano_cavgs_mac/holdout_scoring/<trajectory>/`.

## What was built

- `quality_context=nano` in `model_cavgs_rejection` (UI choice, per-mode check, commander):
  - Gates: empty class, invalid pixels, both local variances at `EPS`, no measurable spectrum.
  - Bank: `lattice_peak`, `spot_count`, `spot_arc`, `spot_diam`, `spot_snr`, `band_presence`.
  - The floor and the guard above; the policy `lattice`.
  - The relational term on a nano band (15 A to 1.5 A, 3 px), because the model schema requires it.
  - `analyze` without a nano model is the hard-gate-only decision. Its table header records the floor and the
    median particle diameter.
  - `apply` with a nano model writes `diam` on each scored class and its median on the class-average entry of
    `os_out`.
- The bank is now a property of the context. Every bank array is indexed by position in the context's own list
  (`cavg_quality_bank_size`, names and families by context). The largest bank, `CAVG_QUALITY_NFEATS`, stays 14,
  because the nano list has six. Extraction, normalization, the table writers and reader, the policies, the model
  file, the scoring and the promotion snippet follow the context. The cryo-EM code paths are unchanged.
- `src/main/cavg_quality/simple_cavg_quality_lattice.f90`, a class (`new`, `extract`, `fit_floor`, `kill`) plus the
  floor normalization. It is ported from the three modules of `src/main/frame_quality/`, which are deleted.
  Dropped as dead: the band maximum, the running median, the white-noise reference and the per-frame table. On
  real class averages its values equal the numpy reference of the previous run within 1e-4.
- The sub-suite "cavg quality nano" of `unit_numerics`, `simple_cavg_quality_lattice_tester.f90`, with 63 checks:
  - lattice against noise; the floor and its guard on mixed, all-lattice and all-empty projects; the gates;
  - the diameter on a 25 A sphere (25.0 A) and disc (28.5 A); `spot_arc` on straight against bent fringes;
    `band_presence` unmasked;
  - the bank by context; the model-file round trips (nano, preset-only); the table round trip through `learn`
    for both banks;
  - the parallel extraction with a team against one thread.

## Defects

Fixed in the class-average program:

1. Project-backed `evaluate` failed for every relational model. It passed no relational feature to the scorer.
2. A model file that names a `preset` without restating its interaction terms was refused, contrary to the
   policy note.
3. `evaluate` accepted training files and a model of different banks. It now refuses them.

Also fixed, a mistake of this run caught by the new tester before any use: a half-plane mask in the ported
lattice code had lost its lower array bound.

Found outside the file table, not fixed:

1. `simple_strategy2D_srch.f90` (`prep4srch`, around lines 154-185). On a fresh 2D start the search sets each
   frame's previous class to a random populated class. The time-chunk assignment written by `init_tseries_refs`
   is therefore never read, and `strategy2D_tseries` moves each frame only among the neighbours of a random
   class. The classes of `analysis2D_nano` are not time chunks: 0.3 to 0.5% of consecutive frames share a class
   and a class spans 0.85 to 0.9 of the trajectory. That holds on the delivered runs, on the previous run's
   `refine2D_nano`, and on this build's `analysis2D_nano` of a hold-out.
2. `simple_nano_density.f90`, line 373 (`profile`), an array bound mismatch (`prof` 448 against 1001) in
   `test_transfer_kernels`. It is the one failure of the Debug fast gate (`unit_single`). The Release fast gate
   passes.
3. The comment on `quality_context` in `simple_parameters.f90` still lists three contexts.

## The three copied files

The commit the worktree was made from (8e0ea79a9) uses three modules that exist only untracked in the maintainer's
checkout. They were copied in verbatim, are not part of this work, were never edited, and are left out of
`full.diff`:

| file | sha256 |
| --- | --- |
| `src/main/nano/simple_atomic_species.f90` | `ee0f0fb947f2f94bbd76bd6577340a1b8ad38217474ce92c47ae3d297160ebbd` |
| `src/main/nano/simple_nano_density.f90` | `8cfdff817664226e9024a64bef8cc5ed47221fff75909cf1027c79ee151572ae` |
| `src/main/nano/simple_nano_density_tester.f90` | `84efcce09666cacf70832deb6d20df2c55a77277148d5819467d829a95bd9596` |

## What was built and run, with results

- Builds:
  - Release, out of tree: gfortran 16, OpenMP, ARCHOPT, OpenBLAS 0.3.27. Its fast gate passes 15 of 15.
  - Debug, bounds-checked: its fast gate passes 14 of 15, the one failure above, before and after the change.
  - `unit_numerics` passes in both builds.
  - The new sub-suite passes, and passes three times at 11 threads, bit-identical to one thread.
  - `check_test_registry.py` passes.
- The cryo-EM byte gate. `learn`, `promote` and `evaluate` run on five synthetic training tables (the repository
  has no class-average fixtures). `analyze` runs in the chunk, pool and sieve contexts on a cryo-EM project with
  a manual selection and on a nano project. All 49 output files are identical before and after the change. The
  newly working project-backed `evaluate` has a baseline of its own.
- The three labeled datasets were rebuilt into projects from the delivered fields (import, `make_cavgs`, class
  states). Populations, states, alignments, `indstk` and `eo` equal the delivered ones class by class and frame
  by frame.
- `analysis2D_nano` ran on all seven hold-outs. `analyze quality_context=nano` ran on the three labeled projects.
  Leave-one-out `learn` and project-backed `evaluate` ran on each labeled project. `apply` ran on copies of the
  hold-outs, and on a copy of a labeled project with `prune=yes`. The prune kept the selected frames of the
  accepted classes, with `indstk` and order intact.
- The by-hand gate was run end to end on a copy of a labeled project, as the policy note writes it: `analyze`,
  the table filter, two `selection` calls. It kept 7,244 frames, exactly the frames of the 378 classes at or
  above 2.09 floor sd.
- Phase 4 (3D): cancelled by ruling. Two of its five refinements had started and were stopped; none finished. No
  3D table exists.

## What is left

- A balanced-accuracy objective for the nano learner, or simply the floor threshold as the product's decision.
  The learner's false-negative-tolerant objective was built for cryo-EM chunks. On nano it picks thresholds too
  low to beat the floor.
- The floor guard on wide floors: two hold-outs with a clear lattice population fail it at 2.6 and 2.9 sd. A
  floor from a three-component fit, or a guard on the gap alone, would deserve a test on more labeled data.
- More labels. Three labeled trajectories are the minimum for leave-one-out; the hold-out tables are ready for
  review.
- The `tseries` fresh-start defect in `simple_strategy2D_srch.f90`, and the bounds error in
  `simple_nano_density.f90`.
- If the class diameter is wanted downstream, `autorefine3D_nano` can read `diam` from `os_out` in place of
  `estimate_diam`. That is a change in the nano3D commander, which was outside this run.
