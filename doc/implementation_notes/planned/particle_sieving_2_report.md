# particle_sieving_2: report

Run of [the plan](particle_sieving_2.md), 2026-11-23, interactive, on the development machine (Ubuntu, gfortran
15.2, 16 cores) with Phase 4 on thunder (gcc-toolset-15, 72 cores); decisions 23 to 26 on 2026-11-23 and 24. Base
`638d99ae3`. Everything is uncommitted for review.

## 1. What was learned

- **A greedy round cannot start fractionally from `solve2D`'s output.** `solve2D` runs its stages at a crop and
  its final `make_cavgs` overwrites the last averages and the FRCs with full-box ones, while the carried class
  sums stay at the crop (MotAB: 320 px averages, 160 px sums). `refine2D` pads smaller references but stops on
  larger ones and treats carried sums of another box as missing. Ruling 15 runs the rounds at the full box: round 1
  is a full update and publishes full-box sums, the later rounds blend them (miniature: fraction 0.706 then 0.300,
  `CLASS SUMS blended ... mass ratio 1.002`).
- **Unseeded, the particles ingested by batches never get a valid sigma2 record.** A pass rewrites a particle's
  record over its band only; the batch particles were inactive when the seed's `calc_pspec` filled every shell, so
  every round's commit skipped them (`canonical sigma2 reduction skipped 600 / 1200 / 1500 active particle
  records`) and the run's noise spectra were the seed particles'; a round without the seed's particles stopped on
  an empty half (a mutation test). By Joe's decision 21 each batch's records are now seeded from the store's global
  curves with `sigma2_state_extend` after its prediction: the warnings are gone, the last round's store holds a
  valid record for every particle (a new gate, which fails with the seeding removed), and the one-thread runs stay
  byte-identical. The extension forms each half's mean over the seeded rows, so a set of one half (halves follow
  micrographs) waits for the next batch. The result still carries no store (decision 20), as
  `particle_sieving`'s has none.
- **Left to the seed's `solve2D`, every batch particle was in the even half** (decision 22). The even/odd
  transaction gives each unmarked stack a half by micrograph group, to the half with fewer active particles; with
  only the seed's particles active, the stacks without active particles add none, so all of them went to the
  lighter half and were marked for good (real set: 50 601 batch particles even, 0 odd). Every round then summed
  them into the even class sums, and no batch could be seeded with sigma2 records (decision 21's first real run
  waited on one half throughout). The program now makes the halves of every selected particle before the seed
  (real set: 233 micrograph groups, 26 569 / 26 481). The miniature hid it (index labels, no micrograph rows); its
  fixture now has a micrograph per stack and no halves, and fails two new gates without the fix. The two Phase 4
  runs of the defaults (section 5) and the first seeding run had the defect.
- **A step's copy must keep the work project's name.** The sigma2 store's layout digest includes `projname`: a copy
  named after its step would have the store refused and rebuilt by `calc_pspec` at every step.
- **A test program named like the program it runs hides it.** In `simple_test_exec`, `params%new` of the program
  found the test's descriptor (no project required) and skipped the project copy. The test is
  `particle_sieving_2_mini`. Also, `simple_test_exec` makes no execution directory for a program it runs in
  process; the program runs in the test's folder.
- **The CNN's threshold must be double precision.** In single precision 0.9 is 0.89999998, so `p_junk` = 0.9 was
  junk against decision 3's "strictly above"; the unit suite caught it on its first run.
- **The CNN fixture cannot exercise the labels at 0.9:** its random weights give `p_junk` 0.48 to 0.52 (MotAB, and
  0.49 to 0.52 on the miniature); the trained `cavg_cnn_v3.bin` puts 86 of MotAB's 100 classes above 0.9.
- **The miniature is fast and reproducible:** 29 to 130 s at eight threads, 92 to 141 s at one, byte-identical at
  one thread over every declared output with no exclusion.
- **On a real set the program runs and keeps its contract.** With decisions 21 and 22 (`run4_mics10`, a
  10-micrograph seed, 10 classes): six batches of about 10 000 and a tail of 656, each with both halves and
  seeded, no sigma2 record skipped, a full round 1 (92 s) and a fractional round 2 (0.025, blended with the
  carried sums, mass ratio 1.009); 32 412 of 53 050 accepted, 0 mismatches with the last labels; 18:45 in all.
  The defaults with both decisions (`run5_defaults`, 24 threads): 4 509 s, seed 3 733 s, student 401 s, the one
  batch (7 089, both halves) 15 s, round 1 357 s; 45 359 accepted, 0 mismatches; junk classes 18 after the seed,
  25 after the round. Before the decisions (the two runs of the defaults, every batch particle in the even half;
  thunder, the 233-micrograph `mic_series_rbp` source, 53 050 particles, box 224, `cavg_cnn_v3.bin`): 5 355 and 5 800 s, of which the seed's `solve2D` 4 595 and
  4 775 s, the student 359 and 674 s, the one batch 17 and 14 s (prediction 0.67 ms and alignment 0.9 to 1.2 ms per
  particle), the full-update round at 224 px 381 and 334 s. 45 580 and 45 411 of 53 050 accepted; junk classes 8
  and 14 after the seed, 17 and 27 after the round. The accepted sets share 40 088 particles (Jaccard 0.79): the
  multi-threaded seed differs between runs, and the labels follow it.
- **The seed dominates and was slow here:** about 2.5 min per `solve2D` iteration against about 33 s in the study's
  run of the same 200 micrographs, which ran alone, at 16 threads, with a mask from its class averages (145 A) where
  this program uses the box's (283 A, decision 2). Not separated further.
- **Under the box's mask the seed's `solve2D` turns many classes into rings** (decisions 23 and 24, 2026-11-23). On
  the betagal stream set its stage 2 (the first centring stage) collapses classes into rotationally symmetric ring
  averages; under the box's mask (283 A), or in the automatic mode, whose stage 2 runs under the box's mask, they
  persist; under 174 to 180 A they recover. `mskdiam` is now taken from the command line, and `mskdiam=0` runs the
  automatic mode. Share of the accepted particles in protein-view classes (judged by eye from the final montages;
  `p_junk` > 0.9): box 41.5 % of 45 359 (run 5), 174.72 A 50.0 % of 36 720 (run 6), 174 A 48.8 % of 36 614 (run 7),
  automatic (140 A) 37.9 % of 50 044 (run 8); `solve2D` alone on all particles at 174.72 A, 46.6 % of 37 688. No
  run rejected a protein-view particle.
- **A 20-micrograph seed is thin, and retraining the student buys nothing measurable** (decisions 25 and 26,
  2026-11-24). Run 9 (`mskdiam=180`, 20 micrographs, 4 177 seed particles, batches of 1 000 by a variant build, a
  round every 5 batches, the student retrained after each round): 3 918 s, of which training 3 424 s; its first
  round relabelled 29 classes. Runs 10 and 11 seeded by the new rule (45 micrographs, 10 083 particles), trained
  once (1 109 s, training 504 s) and after every round (4 008 s, training 3 398 s). The share of new particles a
  round leaves in the student's class: 45 to 65 % once, 45 to 63 % retrained, on the same batches (0.3 points apart
  on average); a replay of run 9's batches with its first student lost 1.7 to 5.5 points against the retrained
  ones. The next round searches every new particle against all classes, so the student only places it until then.
  Accepted (`p_junk` > 0.9; protein-view share): 45 062 (50.0 %), 42 918 (47.6 %), 38 266 (53.5 %) in runs 9, 10
  and 11; the accepted sets' Jaccard 0.71 to 0.75 between any two, as for independent runs.
- **The CNN's most probable group selects more cleanly than `p_junk` > 0.9:** protein-view share 62.0, 63.4, 63.3
  and 68.8 % in runs 7, 9, 10 and 11 against 48.8, 50.0, 47.6 and 53.5 %, still rejecting no protein-view particle.
  The 0.9 rule keeps blurred blobs, edges and noise discs that the CNN calls most likely junk at a lower `p_junk`;
  both rules keep comets and some edges. Decision 3 stands; section 6.
- **The selection moves between rounds by the labels, not the particles.** In run 9, from round 4 on 0.8 to 3 % of
  the earlier particles change class per round (8 to 11 % of those a round searches), while 4 to 11 % change between
  accepted and rejected: 7 to 29 classes change label per labelling, most with `p_junk` between 0.8 and 1.0 on both
  sides of 0.9. The result takes the last round's labels.
- **The first greedy round empties the seed's smallest classes for good:** 17, 17 and 23 of 100 classes in runs 9,
  10 and 11, the seed's empty or nearly empty ones (up to 21 particles); nothing refills a class, so these runs sort
  into 77 to 83 classes. Section 6.
- **A run is reproducible only with `SIMPLE_SEED`:** without it the seed draws from `/dev/urandom`; runs 10 and 11
  put 75 of the same 10 083 seed particles in the same class. The plan's identity is at one thread with a fixed seed.

## 2. What changed, phase by phase

| Phase | Files | What |
| --- | --- | --- |
| 0 | 1 | the plan: decisions 15 to 20, Progress |
| 1 | 1 | `simple_cavg_cnn_predict.f90`: `cavg_cnn_class_probs` extracted, the program's table and log unchanged |
| 2 | 9 | `src/main/sieve/simple_ptcl_sieve2.f90` and `_tester.f90` (decision 21: `seed_sigma2`, `both_halves`), `src/defs/simple_defs_sieve2.f90` (new); `commander_sieve_ptcls2` in `simple_commanders_sieve.f90`; the route in `simple_exec_sieve.f90`; the descriptor in `simple_ui_preproc.f90`; the suite row in `simple_commanders_test_class.f90` and its `suite=` help in `simple_test_ui_class.f90` |
| 3 | 7 | `simple_commanders_test_highlevel_sieve2.f90` (new) and its declaration in `simple_commanders_test_highlevel.f90`; the route in `simple_test_exec_highlevel.f90`; the descriptor in `simple_test_ui_highlevel.f90`; the registry row; the CTest entry and `SIMPLE_CTEST_BUDGET` 46 in `production/CMakeLists.txt`; the budget's number in `test_environment_policy.md` |
| 5 | 3 | section 15 of `ptcl_sieve_policy.md`; this report; the notes index |
| decisions 23, 24 | 5 | `mskdiam` captured by the commander before `params%new` (0 dropped: the seed's estimate); `mask_diameter` and the seed's registered mask in the module; its test; the UI help; the policy |
| decisions 25, 26 | 8 | `SIEVE2_SEED_NPTCLS` and `SIEVE2_NMICS` 100 in `simple_defs_sieve2.f90`; the seed's two limits in `plan_batches` and one training in `run`; `test_plan_seed_cap`; the workflow test's seed ended by a cap of 450 before `nmics` 10; the UI help and default; the policy; the plan; this report |
| decision 40 | 8 | `mrc2jpeg_tiled`'s optional `frames` (a true-colour montage with red and orange frames) in `simple_gui_utils`; `test_framed_tiles` in the JPEG tester; `frame_codes` and the two montages in the module, `test_frame_codes`; the workflow test checks them; the policy; the plan; this report |
| decision 39 | 7 | `SIEVE2_ROUND_CAVGS_DIR`, `SIEVE2_ROUND_CAVGS`; `round_cavgs` (`make_cavgs` on a frozen round's particles) in the module; the workflow test's freeze gate checks the stacks; the registry declares them; the policies; the plan; this report |
| decisions 37, 38 | 8 | one labelling in `run` (at the freeze, else after the last round, or on the seed's averages), the reports without labels (`pct_str`, - columns, the run total over labelled rounds), the log line without empty parts; `restore_cavgs` (default no in the commander, `final_restore` in the settings, the UI input); the workflow test's one-labelling gate and its restoring last round; the policies; the plan; this report |
| decision 36 | 7 | `SIEVE2_FREEZE_CORR`, `SIEVE2_FREEZE_NROUNDS` and their settings; the freeze count in `run`, the frozen round's new-only fraction and `restore_cavgs=no` in `greedy_round` (`freeze_count`, `newonly_frac`); `test_freeze_rules`; the workflow test's freeze gate (frozen after one round); the policies; the plan; this report |
| decision 35 | 6 | the class averages' change (`refs_change`, `class_pops`), the column and log line; its test and gate; the policies; the plan; this report |
| decisions 33, 34 | 9 | the agreement split by the predicted class's label (two columns, the log line, the run total); `student_by_class.txt` (`class_agreement`, `kept_order`); their tests and gate; the table among the identity outputs; the policies; the plan; this report |
| decisions 31, 32 | 8 | `SIEVE2_ROUND_OLD_FRAC` and the settings' `round_old_frac`; `greedy_frac` with the earlier particles' share; the round's agreement counts (`round_agreement`), the log lines, the run report's two columns and the run total; `test_greedy_frac` rewritten, `test_round_agreement`; the workflow test's agreement gate and reads; the two policies; the plan; this report |
| decision 30 | 7 | `SIEVE2_LPSTOP`; `lpstop` defaulted and refused (not positive) by the commander, handed to the seed's `solve2D`, and each round's `lp` from the class FRCs held at it (`round_lowpass`) with `stage`/`gold_std_stage`; `test_round_lowpass`; the workflow test's gate on the rounds' limits; the UI input; the two policies; the plan; this report |
| decisions 27, 28 | 9 | the key `retrain` in `simple_parameters.f90` and `_parse.f90`; the commander (parsed and refused before `params%new`; `mskdiam` required); `retrain_rounds` in the settings, `retrains_after` and the warning for a round no batch follows in the module, the folders `learn_<nn>` again; `test_retrains_after`; the workflow test retrained after round 1 (a training-folder gate); the UI input and `mskdiam` required; the policy; the plan; this report |

## 3. Decisions taken under delegation

- 18, the rounds' iteration numbers (`SIEVE2_ROUND_ITER0` = 8): past `refine2D`'s onsets of the in-plane shift
  search and of the band from the class FRCs, as `solve2D`'s terminal pass is.
- 19, a CTest entry of its own (budget 45 to 46), **for Joe's confirmation**: the plan registers the test, and
  `test=identity` finds registered entries only.
- 20, no sigma2 store in the result (section 1); its reason replaced by decision 21 (Joe), which seeds the batch
  particles' records.
- 22, the halves made by the program over every selected particle before the seed (section 1).
- Without a number: every step's copy under the work project's name; the test's name; the threshold in double
  precision; the test's outputs relative to its own folder.

## 4. Deviations from the plan

- **The children's two base runs at one thread** (Phase 0) were taken inside the program's own one-thread runs, not
  on their own: no child changed in this run but `cavg_cnn_predict`, which was checked byte-identical in Phase 1.
- **The Debug build** was a separate tree (`build_debug`, removed afterwards), not `compile_debug.sh`, which would
  have replaced the Release tree.
- **Default f** ("the seed's last stage's settings") became the greedy round's own command line: the seed's last
  stage runs at a crop and uses `solve2D`'s stage keys; the round sets `refine=greedy`, `autoscale=no`,
  `center=no`, `startit`, `minits=maxits=1` and `update_frac` from round 2, and takes `refine2D`'s defaults for the
  rest (`objfun=euclid`, `ml_reg=yes`, measured sigma2, `trs=MINSHIFT`, the band from the class FRCs).
- **Batches of 300 in the miniature**, not 250: the fixture's stacks hold 100 particles, so a batch closes at the
  first stack that reaches 250, and the expected counts are exact at 300.
- **Files outside the table** (a row added on 2026-11-23): `src/main/exec/simple_test_exec_highlevel.f90`,
  `production/CMakeLists.txt`, `doc/policies/test_environment_policy.md`. **Not needed:** `src/main/ui/simple_ui.f90`
  (its `particle_sieving` line is the stream's JSON), `sieve.inf` (sources are globbed).
- **Runs 9 to 11 ran variant builds on thunder** (`~/Code/SIMPLE_ps2_b1000`, not this tree): batches of 1 000 (the
  constant changed), and for runs 10 and 11 the seed rule and an environment switch for the one training. The
  module now takes run 10's path (the seed rule, one training) without a switch.

## 5. Evidence

| Gate | Result |
| --- | --- |
| base build and fast gate | Release; 15 of 15 |
| Phase 1, `cavg_cnn_predict` table | byte-identical to the base on MotAB's 100 averages, fixture and trained weights |
| `cavg_cnn`, `cavg_cnn_input` | 41 of 41, 12 of 12 (Release and Debug) |
| `unit_project suite=particle_sieve2` | 26 of 26 (Release; 22 of 22 in Debug before decision 21; 36 of 36 with decision 25's test); its first run caught the precision defect |
| `unit_project suite=particle_sieve` (unchanged sieve) | 119 of 119 (Debug) |
| `unit_ui`, `check_cli_chain.py`, `check_test_registry.py`, `check_submodule_linkage.py`, `check_guide_imports.py`, `check_descr.py` | pass |
| fast gate, final code | 15 of 15 (Release); 15 of 15 in Debug (over its 60 s time budget, as Debug is) |
| `particle_sieving_2_mini`, eight threads | 9 of 9 gates in four runs (Release, CTest included) and under Debug bounds checks; with decisions 21 and 22, 12 of 12 (the sigma2 records' and the two halves' gates added) |
| identity, one thread | 14 of 14 declared outputs byte-identical, four times (before and after decisions 20, 21 and 22) |
| mutation | rounds without batch 1's particles and every particle accepted: the round, active-particle and accepted-particle gates and the unit suite fail; the seeding removed: the sigma2 records' gate fails and the rounds skip 600 / 1200 / 1500 records again; the halves left to the seed: the batches' halves and sigma2 gates fail (every batch waits on one half) |
| Phase 4, thunder, two runs (before decision 22) | normal stop; 53 050 decided, 0 mismatches with the last labels, the round's copy 53 050 active, one batch (7 089) and one round |
| thunder, `run4_mics10` (decisions 21 and 22) | normal stop; 0 mismatches; halves 26 569 / 26 481; six batches with both halves, all seeded; no sigma2 record skipped; rounds full and fractional (0.025) |
| thunder, `run5_defaults` (decisions 21 and 22) | normal stop; 0 mismatches; the batch with both halves and seeded; no record skipped; 45 359 of 53 050 accepted |
| thunder, runs 6 to 8 (decisions 23, 24) | normal stops; masks 174.72, 174 and 140 A (automatic); protein-view shares in section 1 |
| thunder, runs 9 to 11 (decisions 25, 26; variant builds) | normal stops; 3 918, 1 109 and 4 008 s; section 1 |
| decisions 25 and 26, this tree | `particle_sieve2` 36 of 36; `particle_sieving_2_mini` 12 of 12 at eight threads (21 s) and twice at one (81 s), the seed ended by its particle cap; identity 14 of 14; `unit_ui`, the checks and the fast gate (15 of 15) pass |
| decision 40, this tree | `mrc2jpeg` 33 of 33 (the frames decoded); `particle_sieve2` 58 of 58; `particle_sieving_2_mini` 17 of 17 at eight threads and twice at one; identity 15 of 15; the fast gate 15 of 15 |
| decision 39, this tree | `particle_sieve2` 57 of 57; `particle_sieving_2_mini` 17 of 17 at eight threads and twice at one (the frozen rounds' stacks); identity 15 of 15; the fast gate 15 of 15 |
| decisions 37 and 38, this tree | `particle_sieve2` 57 of 57; `particle_sieving_2_mini` 17 of 17 at eight threads and twice at one (one labelling after round 1; the last round restoring); identity 9 of 9; `unit_ui` 1 042 of 1 042; the fast gate 15 of 15 |
| decision 36, this tree | `particle_sieve2` 57 of 57; `particle_sieving_2_mini` 17 of 17 at eight threads and twice at one (frozen from round 2); identity 12 of 12; the fast gate 15 of 15 |
| decision 35, this tree | `particle_sieve2` 51 of 51; `particle_sieving_2_mini` 16 of 16 at eight threads and twice at one; identity 15 of 15; the fast gate 15 of 15 |
| decisions 33 and 34, this tree | `particle_sieve2` 50 of 50; `particle_sieving_2_mini` 16 of 16 at eight threads and twice at one; identity 15 of 15; the fast gate 15 of 15 |
| decisions 31 and 32, this tree | `particle_sieve2` 44 of 44; `particle_sieving_2_mini` 15 of 15 at eight threads and twice at one, the agreement reported for every round; identity 14 of 14; `check_cli_chain.py` and the fast gate (15 of 15) pass; `retrain` with a single round fixed and checked by hand |
| thunder, runs 12 and 13 (this tree, batches of 1 000) | normal stops; 510 and 544 s; every round at 10 A; the student 58 % in all; retraining after round 5 about +3 points in the next two rounds |
| thunder, run 14 (a variant: a greedy round on the seed before the training) | normal stop, 502 s; the student 57.4 % against 58.1 %: not taken |
| thunder, runs 15 and 16 (variants: the freeze, then with new-only rounds) | normal stops, 520 and 483 s; frozen from round 5; 16 and 5 labels changed at the end; run 16's behaviour taken as decision 36 |
| decision 30, this tree | `particle_sieve2` 38 of 38; `particle_sieving_2_mini` 14 of 14 at eight threads and twice at one, the rounds at 10 A (6.4 to 5.7 A from the FRCs before), the seed's e/o stage at 6.4 A as ruled; identity 14 of 14; `unit_ui`, `check_cli_chain.py` and the fast gate (15 of 15) pass |
| decisions 27 and 28, this tree | `particle_sieve2` 39 of 39; `particle_sieving_2_mini` 13 of 13 at eight threads (27 s) and twice at one (82 s), retrained after round 1, the last round of the list ignored; identity 14 of 14; `unit_ui`, `unit_core suite=parameters` 66 of 66, `check_cli_chain.py` and the fast gate (15 of 15) pass; the refusals of a missing `mskdiam` and a bad `retrain` checked by hand |

## 6. Left open for the maintainer

- **The sigma2 seeding's one-half wait** (decision 21) is covered by the unit suite only: every batch of the
  miniature (three micrographs) and of the real runs held both halves.
- **Reported, not gated:** the real set's results on thunder, `~/ps2_phase4/run5_defaults/report/` and
  `run4_mics10/report/` (the last round's montage, averages and halves, the labels tables, the run report, the
  particle fields; `run1`, `run2` and `run3_eo_one_half` are the runs before decision 22), and whether the CNN's
  `p_junk` at 0.9 is a usable rejection criterion: the trained model marks 86 of MotAB's 100 classes junk, and on
  the thunder set 18 of 100 after the seed and 25 after the round (every run's count grows after the round).
- **Not exercised:** `nparts` (the children distributed). A fractional round on real data was exercised only with
  a 10-micrograph seed (`run4_mics10`); the defaults give this set one batch.
- **The seed's time on thunder** (about five times the study's per iteration): whether the box's mask, the two
  concurrent runs or the load caused it is not separated.
- **The key `retrain`** (decision 27), a new parameter, was approved by Joe on 2026-11-24.
- **Ruled 2026-11-24 (decisions 28 and 29):** `mskdiam` required; `p_junk` > 0.9 kept although the most probable
  group was cleaner on every run read (section 1); classes the first round empties stay empty; every step keeps its
  copies (28 and 24 MB each at 53 050 particles; with batches of 1 000 a run took 4.2 to 6.3 GB); batches of
  10 000 and a round every 5 kept, which gives betagal one round, after the last batch; decision 19 confirmed.
- **Addition to the `simple-main-sieve` guide** (applied 2026-11-24 on Joe's word): "`particle_sieving_2` (`simple_ptcl_sieve2`)
  is a separate controller: a seed `solve2D`, the class-average CNN's labels, the class student trained on the seed
  (again after the rounds `retrain` lists) and used for the batches, and a greedy `refine2D` round at the full box
  every five batches. Each step runs on a copy of the work project under
  the project's own file name in its own folder (the sigma2 store's layout digest holds `projname`); the result
  carries no sigma2 store. Policy section 15; tests `unit_project suite=particle_sieve2`, CTest
  `particle_sieving_2_mini`."

## 7. Commit message

```
particle_sieving_2: a seed solve2D, CNN junk labels and the class student for the rest

A new program beside particle_sieving (simple_ptcl_sieve2, a controller of whole
programs): solve2D on the particles of the first micrographs, up to nmics (100) or
10 000 particles, the class-average CNN's labels (junk above p_junk 0.9), the class
student trained on the seed, the rest predicted and aligned in batches of 10 000,
and a greedy refine2D round at the full box every five batches (the first a full
update, the later ones blended with the carried class sums), each followed by the
labels and, after the rounds the new key retrain lists, the student's retraining.
Particles of junk classes are deselected at the end; the result carries no sigma2
store. mskdiam is required and taken by every step (0: the seed's solve2D estimates
it). lpstop (10 A by default) holds the alignment: the seed's solve2D marches to it
and each round's limit from the class FRCs is held at it. The program makes the
even/odd halves of every selected particle
before the seed (left to the seed's solve2D, every batch particle fell into one half)
and seeds the batch particles' sigma2 records from the store's global curves
(sigma2_state_extend). cavg_cnn_class_probs extracted from cavg_cnn_predict (its table
byte-identical). Unit suite particle_sieve2; CTest entry particle_sieving_2_mini
(SIMPLE_CTEST_BUDGET 46), byte-identical at one thread.
```
