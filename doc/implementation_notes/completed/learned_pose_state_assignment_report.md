# Learned pose assignment: report of the run `pose_student`

Report, 2026-11-08. The run carried out the plan
`doc/implementation_notes/completed/learned_pose_state_assignment.md` on the dell: the home build
machine, 24 cores of an Intel Core Ultra 9 285, 62 GB, no GPU, Oracle Linux 8.10 and gfortran 15.2.
It ran unattended from 2026-11-06 to 2026-11-08 as one uncommitted diff against HEAD d4e89630c (the
multi-state posterior draw). Each phase's diff is `R/review/phase_N.diff`, written by the run's
driver when the phase ended; `R` below is `/home/elmlundho/agent_runs/pose_student`. The plan's
Progress section holds the detail and every evidence path; this report is the summary.

**Changed after the run (2026-11-09):** the student now classifies the 500 mirror pairs of the grid
and HHCNN is invariant to the mirror, the search resolving the pair through the Euler relation; on
the medium assembly's whole training half it puts 94.0% of the centred test half within 6° of the
truth against 91.5% for the 1,000-class student described below. The plan's Progress entry "Mirror
pairs (2026-11-09)" holds the change and its measurements; this report describes the run as it ended.

The change builds the pose student: a small network, HHCNN, the Harmonic Hermitian convolutional
neural network (Hermitian because every layer commutes with conjugation on the half ring with the
conjugate wrap, so the network is a model of the Fourier transform of a real image; harmonic because
the readout's class templates live in the angular harmonic domain of the polar transform and the
rotation search is the cross-correlation over the cyclic shifts through the angular Fourier
transform; convolutional on the polar Fourier representation), trained per dataset on the
probabilistic search's own posterior rows over a grid of 1,000 directions, that proposes for a
particle the search has not seen a direction, a cone about it and a confidence, and hands them to
`refine3D`'s neighbourhood search (`refine=prob_neigh prob_neigh_mode=geom`), which finds the
in-plane rotation, the shift and the pose with its own machinery.

- Two private programs: `learn_pose_student` (the teacher's rows from the production fill at the
  student's grid and band, the training of the designs cheapest first, the model file) and
  `predict_poses` (the per-particle direction, cone and route into three new orientation slots).
- The neighbourhood search takes a particle's cone from its record, about the stored direction and
  about its mirror partner; a particle without the slot runs the base's code path.
- Validation: the split-half protocol on the Ribosembly truth sets, vcrnapol and seven real sets;
  the training-set size curves; the cost table; a nightly route (`pose_student`); the algorithm page
  `doc/algorithms/pose_student.md`.

No parameter was added: the band is `lp`, the grid `nspace`, the model `infile` and `fname`, the
cone a numbered orientation slot. The band and the width are the only hyperparameters.

## What we learned

- **The existing routes are untouched, and the machinery works.** A project without student
  records runs the base's code: the `simulated_workflow_1jxy` fixture from the base build and from
  the tree gives 216 declared files the same, and three `prob_neigh` iterations of trpv1 (35,645
  particles, four parts) give every volume, half map, FSC file, reprojection model and sigma2 state
  byte for byte, the particle fields identical by their printouts. The fast gate passes on every
  phase's tree (15 of 15), the three new unit suites pass (polar features 39 checks, HHCNN 26, pose
  student 34), and the nightly tier on the final tree passes every one of its 20 entries.
- **The teacher is not the cost.** The dense fill at the student's grid in the 15 Å band costs
  1.9 to 2.2 ms per particle wall on 24 threads: the training half's rows are a minute. The fill
  the plan feared would be slow is not the problem of this design.
- **The student's input and output are right, its first head was not.** The plan's invariant head
  (the angular mean of the real part and of the modulus) sees 8° of the ring before pooling and
  learns the marginal over directions; extended to the angular harmonic magnitudes it reached a
  held-out loss of 4.39 (uniform 6.91) and 10.6% exact mode agreement on the medium assembly, and
  its cones held the truth for 85% of centred unseen particles. Those cones were the whole sphere:
  the posterior's 99% set held the mode's neighbourhood and its mirror partner (the antipodal
  direction, whose projection is the reflected image), and the seeded rows' better reads were an
  exhaustive search over 20,000 directions, not a localisation. Hans's ruling R1 folded the
  posterior over mirror pairs and bounded the route to 10°; under it that student seeded nobody.
  Ruling R2 found the student is not shift-invariant (the real-part mean and the harmonic
  magnitudes carry the shift; a modulus-only head is invariant and learns nothing). Ruling R3
  replaced the head with a rotation-searching template readout: one chiral template per direction
  class in the angular harmonic domain, the score the maximum over the ring's rotation grid. With
  it the medium assembly's student reaches a held-out loss of 1.78 and a mode agreement of 46%
  (5% before), recovers 45% of the held-out seeds under the folded rule (0% before), and seeds 46%
  of the centred test half with cones of 6.5° to 9.5° that the neighbourhood search aligns as well
  as the dense search aligns the rest (99.0% against 98.7% within 6° of the truth after `cont`) at
  a third of its cost per particle (3.2 against 9.8 ms), the prediction costing 0.56 ms per particle
  against the coarse fill's 2.2 ms. That is the student the plan proposed, on the set the plan
  proposed it for.
- **On the real sets at 15 Å the student seeds almost nobody, because its teacher is broad.** One
  dense iteration at 15 Å from nothing puts 6 to 46% of a real set's particles within 6° of the
  set's reference alignment (trpv1 34%, bgal 24%, streptavidin 6.5%, CESPED TEST 36.5%, FlhB 29%,
  11120 46%; 10648 69%, the one sharp teacher). The student learns that breadth: its held-out mode
  agreement rises from 3-5% with the invariant head to 10-25% with the readout, but no held-out
  particle gets a confident 99% set within 10° of its mode, so the threshold admits nobody, and the
  protocol costs the dense search plus the training. On 10648 (PKM2, D2) the student seeds 2 to
  4.4% of the test half from 180 training rows on, those rows aligned as the dense search aligns
  the rest. The student can only hand over what its teacher resolves; the medium assembly's
  teacher, the same search against a map from truth poses, is sharp enough.
- **How few particles.** On the medium assembly with the readout a tenth of the training half
  (4,500 rows) recovers more held-out seeds than the whole half (63% against 45%, with 64 channels
  kept by the rule) and seeds 1.7% of the uncentred test half at 90% within 6° against the dense
  search's 66%: the training budget's early stopping (the best epoch is the third of six on the
  whole half) is tuned for small sets. On the real sets the curves at 200, 1,000 and 5,000 rows seed
  nobody (trpv1) or a few percent at every size (10648).
- **The student needs centred particles.** A student trained on centred particles tolerates about
  a pixel of shift at 15 Å; the 2 px augmentation disk buys a pixel at no cost, a 10 px disk (the
  benchmark's picking error) defeats the fit. On the uncentred medium test half the readout's
  student seeds 1.0% (at 98.7% within 6°), on the centred one 46%. The shift is the open question of
  the design; the search's own shift posterior is the centring record, as ruled.
- **Reconstructing from the student's directions alone** (ruling R4: every row at the student's
  folded mode, the cone search for the in-plane rotation and the shift, one `cont` iteration) gives
  maps as good as or better than one dense search's wherever the particles are centred or the set is
  large: trpv1 at the reference shifts 0.143 at 5.4 Å against 6.5 Å, CESPED TEST at the shifts 9.2
  against 10.3 Å, FlhB 5.9 Å both, 10648 from a 900-row student 3.82 against 3.96 Å, 11120 4.39 Å
  both; worse on uncentred particles (the medium assembly blank 40.8% within 6° against 66.9%). The
  read decides nothing by itself; the maps are the maintainer's to judge.
- **Costs.** The prediction's features and forward pass cost 0.5 to 0.9 ms per particle. The cone
  search on narrow cones costs a third of the dense search per particle; on the wide cones the real
  sets' students give it is a near-exhaustive search over the 20,000-direction grid and costs 2 to 4
  times the dense search (trpv1 12 minutes against 3.4). The training of the designs is the cost of
  the cycle: on trpv1 16 minutes for the whole half, 2.5 minutes at 900 rows, against 173 to 205 s
  for one dense iteration on the test half.

## What changed, phase by phase

**Phase 0, foundations and the teacher's rows** (`R/review/phase_0.diff`, 2026-11-06). The folder
`src/main/learn/` with `simple_polar_feats` (the input tensor from the polar calculator through two
new public getters of `polarft_calc`, the conjugate wrap, the rotation by cyclic shift, the shift
through the calculator's phase ramps, a feature cache under the plane cache's contract) and its
suite; the dense table's row getter `get_posterior_row` and the row record in the table utilities;
the strategy `simple_strategy_pose_student` with the private program `learn_pose_student` in its
rows-only form; the medium assembly's view-balanced halves and the teacher's rows at 15, 10 and 8 Å
(mode within 6° of the truth for 99.3%, the truth inside the 64-entry cone for 100%); the
deterministic in-plane alternative measured (88.6% against 78.1% exact mode, all neighbour swaps);
CESPED TEST imported, 10648 and 11120 downloaded.

**Phase 1, the network and the training** (`phase_1.diff`, 2026-11-06). `simple_circ_cnn` (now
`simple_hhcnn`) with its suite (the Hermitian contract per layer, the CTF normalisation, the gradient
check, the thread identity at 4 threads bit for bit, the round trip, the toy problem at 0.988);
`simple_pose_student` (the model, the training with every tenth row held out, the cone, the file);
`predict_poses` in its first form; the first size-curve points at 15 Å; the shift read. Three
defects found on the real rows: the first block's normaliser took the shell weights twice; the
invariant head is a bag of local statistics (the harmonic magnitudes were added, a departure from
section 6 recorded for the ruling); the model header's six digits.

**Phase 2, the programs and the hand-off** (`phase_2.diff`, 2026-11-07). The three orientation slots
`pconf`, `cone`, `route` (57 to 59 of the record's 64 reals, so the project file's width and every
earlier project read unchanged); the neighbourhood search's per-particle cone; the designs cheapest
first with the resampled gain rule and the confidence threshold on the held-out rows; `predict_poses`
complete (the route, a part's document, the project copy's segment); the identity against the base
export; the protocol on the medium assembly at 15, 10 and 8 Å, the size curve, the small assembly,
vcrnapol; the finding that every cone spanned the sphere.

**Phase 3, real data and cost** (`phase_3.diff`, 2026-11-07 to 08). Rulings R1 to R4 applied: the
fold and the route bound, the shift check in `predict_poses`, the template readout (model files
version 2), the student-only alignment script; the protocol on trpv1, bgal, streptavidin, CESPED
TEST, FlhB, 10648 and 11120 with the size curves on trpv1 and 10648; the cost table from three runs.

**Phase 4, the nightly route and the documents** (`phase_4.diff`, 2026-11-08). The route
`pose_student` (the protocol in miniature on 2,000 blob-phantom particles with their truth poses:
eight machinery checks, the reads reported with no floor by ruling R3; 6:49 at 8 threads), its
registry row, CTest entry, test UI help and router case, the policy's card and budget (43);
`doc/algorithms/pose_student.md`; the plan's sections brought to what was built; the proposal on the
in-loop consumer; the nightly tier green.

**Phase 5, this report** (`phase_5.diff`, 2026-11-08). Ruling R5: the network named HHCNN, the module
and its tester renamed `simple_hhcnn` and `simple_hhcnn_tester` (the type `hhcnn`, the suite
`HHCNN`, the suite key in the test UI), the name defined in the algorithm page, the plan and here;
the fast gate and the three suites rerun; the plan moved to `completed/`.

## Decisions taken under delegation

Each is recorded in the plan's Progress entry of its phase unless another place is named.

- The balanced selection needs a score on every particle and `refine=eval` writes none, so a
  `refine=greedy_inpl maxits=1 trs=0` pass supplies it; the truth poses are then restored from the
  printouts, because `refine3D`'s reference centring (on by default) had moved every stored shift.
  `learn_pose_student` runs its references with `center=no` (phase 0; section 8).
- The row record keeps at most 256 entries with the tail's mass (no row needed more); the feature
  cache follows the plane cache's contract (phase 0).
- The private key dictionary lacked `pgrp`, `fname` and `part`; they were added as existing global
  keys, not new ones (phases 0 and 2).
- The first block is normalised by the ring-averaged CTF squared alone, the tensor carrying the
  shell weights already (phase 1; section 6 as built).
- The training budget: 2,500 optimiser steps in batches of a sixteenth of the rows within 16 and
  256, patience 250, the best epoch kept, so every point of the size curve gets the same steps
  (phase 1). The augmentation disk is three times `shift_sigma`; the measurements used 2 px, the
  plan's 10 px disk learning nothing (phase 1; section 5 as built).
- The route rule: a particle is seeded when its folded 99% set is within 64 directions, its radius
  within the bound (`prob_athres` or the grid's spacing, whichever is larger: 10° at the defaults)
  and its confidence at or above the threshold; the cone's radius is at least the grid's spacing
  (phase 2, ruling R1; sections 7 and 8).
- The width rule's resampling (1,000 draws by micrograph with a Lehmer generator from the run seed,
  the wider design kept when the 5% quantile of the mean gain is above zero) and the threshold (the
  longest confidence-sorted prefix of the held-out rows with 95% recovered) (phase 2; section 7).
- The readout's constants: 8 read channels, 8 ring bands, angular order 16, the templates
  initialised from the classes' reprojections with the scale set so a class's own reprojection
  scores logit 8 (phase 3, ruling R3; section 6).
- The shift check reads the first 256 selected particles at 2, 5 and 10 px in four directions
  (phase 3, ruling R2; section 4).
- The real sets' preparation: the RELION 3.0 STAR files rewritten in the 3.1 layout for the import;
  FlhB's per-micrograph stack names mapped onto its consolidated stack, the maintainer's defocus
  table confirming the order; FlhB and 11120 trained on a class-balanced 20,000 of the training
  half; the CESPED entries keep the benchmark's own halves (phase 3).
- The nightly route runs on the blob phantom of the multi-state fixture, not on the 1JYX set the
  plan named: a D2 tetramer at `c1` gives every posterior four equivalent modes, so the student seeds
  nobody and the seeded-row checks never run (phase 4; section 10 tier 6). The route sets
  `shift_sigma` 2 (the default 10 px is a 23 px disk on a box of 64) and resets the classic
  generator to the gate's seed before each in-process search. Its reads carry no floor (ruling R3,
  item 7); the policy's card says a 20-seed profile sets them.
- The in-loop consumer (section 9): not warranted by the measured cycle as the programs stand; a
  candidate only where the particles are centred by an earlier stage, the teacher is sharp and the
  training is held to the cheapest design at a few thousand rows. The decision is the maintainer's;
  the sampling and fractional-update policy is unchanged (phase 4).

## Deviations from the plan

- The invariant head of section 6 was replaced twice: by the harmonic magnitudes in phase 1 (the
  plan's head learned the marginal) and by ruling R3's template readout in phase 3 (the harmonic
  head's cones spanned the sphere). Section 6 now describes the readout.
- The cone is read on the folded sphere (ruling R1), not on the plan's whole sphere; the route is
  bounded by a radius, not by the count alone.
- The augmentation disk is 2 px in every measurement, not the plan's three sigma of the search's
  centring width.
- The nightly route uses the blob phantom (above), and reports its reads instead of gating them.
- The plan's stricter variant of the protocol (a de novo `solve3D` of the training half as the
  references) was not run: the rulings' reruns took the time.
- The identity against the base export was run in phases 2 and 3 (the folded-cone tree) and not
  again after the readout, as ruling R3 allowed: the readout touches no existing route.
- The student-only table of ruling R4 reports the direction agreement within 6° only; the ruling
  asked for 6° and 12° on the truth sets. The 12° reads are in the runs' tables under
  `R/scratch/keep/<set>/protocol/so_*`.
- The parameters tester, listed for phase 4, was not edited: it lists no test names.

Base behaviours met and recorded, not changed: `refine3D refine=eval` writes no score; the
reference centring of `refine3D` moves a truth set's stored shifts; `import_starproject` resolves no
bare stack name and mishandles `import_dir`; `fname` must carry a known extension; the base build's
project writer, knowing 56 slots, drops the student's three (every project-writing step of a chain
uses the tree's binaries); two 24-thread jobs side by side oversubscribe the machine's spinning
OpenMP teams by 7 to 100 times.

## Evidence against the baseline

| check | where | result |
| --- | --- | --- |
| `simulated_workflow_1jxy` from the base build and the tree, one thread, fixed seed | `R/scratch/identity/1jxy_identity_early.log`, `1jxy_identity_p3.log` | 216 declared files the same, 0 differing (phases 2 and 3) |
| three `prob_neigh` iterations of trpv1, four parts, `SIMPLE_SEED` fixed | `R/scratch/identity/trpv1_identity.txt` | every volume, half map, FSC file, model and sigma2 state byte for byte; the particle fields identical |
| the fast gate on each phase's final tree | `R/scratch/logs/build_final.log`, `build_gate_final.log`, `fast_gate_phase2.log`, `build_p3_readout.log`, `build_p4_final.log`, `build_p5c.log` | 15 of 15 every time |
| the learn suites alone, 16 threads | `R/scratch/unit/p5c/` | polar features 39, HHCNN 26, pose student 34 checks |
| `test=nightly ncunits=24` on the phase 4 tree | `R/scratch/nightly/nightly.log`, `nightly_report.tsv` | 20 of 20 entries passed, no base failure, 46 minutes |
| the `pose_student` route by hand and inside the nightly | `R/scratch/pose_route/route.log`, `metrics.tsv` | PASS; the same reads to the digit |

## Measurements

The split-half protocol: a set's halves are view-balanced; the training half's map is the
references, the student trains on its rows, the test half enters blank (no poses, no shifts) and,
where the halves carry shifts, centred at them; the reads are the final direction's agreement with
the reference alignment (the truth on the truth sets) within 6°, per route, beside the dense
search's on every row. All at 15 Å with the readout unless marked.

| set | training set (a tenth held out) | held-out mode agreement | seeded of the test half | seeded rows within 6° / the dense search's rows | dense search, every row, from blank |
| --- | --- | --- | --- | --- | --- |
| medium assembly, centred | 24,000 | 46.0% | 46.1% | 99.0% / 98.7% | 98.9% |
| medium assembly, blank | 24,000 | 46.0% | 1.0% | 98.7% / 66.6% | 66.9% |
| medium assembly, blank, a tenth of the half | 5,000 | 54.4% (64 ch) | 1.7% | 90.2% / 66.4% | 66.9% |
| small assembly, blank | 12,000 | 20.8% | 0 | - | 47.9% |
| vcrnapol state 1, blank (simulated, SNR 0.1) | 8,500 | 78.2% | 86% | at the dense search's accuracy | 99.99% |
| trpv1, blank and centred | 17,752 | 11.6% | 0 | - | 34.3%, 42.5% |
| bgal | 2,430 | 23.5% | 0 | - | 24.4%, 35.0% |
| streptavidin | 5,701 | 3% | 0 | - | 6.5%, 13.2% |
| CESPED TEST | 5,000 | 10.0% | 0 | - | 36.5%, 42.3% |
| FlhB | 20,000 of 81,204 | 11.6% | 0 | - | 29.4% |
| 10648 (PKM2, D2) | 200 / 1,000 / 5,000 / 20,000 | 0 / 11 / 8.6 / 7.9% | 3.0 / 4.4 / 3.5 / 2.0% | 66.4 / 68.3 / 68.6 / 77.3% against 68.7-68.8% | 68.9% |
| 11120 (TSHR) | 20,000 of 194,019 | 17.7% | 0 | - | 45.8% |
| blob phantom (the nightly route, 8 Å) | 1,000 | 68% | 1.5% | 26.7% / 33.3% | 3.8% |

The cost table (trpv1, the test half of 17,752 particles at box 256, 24 threads, the machine idle;
wall from three runs, core time as user plus system, peak memory as the maximum resident set):

| step | wall | core time | peak memory | per particle |
| --- | --- | --- | --- | --- |
| the rows (the fill at 1,000 directions, the training half) | 39.7 s | - | - | 2.24 ms |
| the training (two designs, early stopping) | 961, 964, 948 s | 19,679 s | 2.7 GB | 60 ms wall per training row |
| the prediction | 17 s | 345 s | 1.6 GB | 0.69 ms (features and forward pass) |
| the dense search on every row | 173, 204, 204 s | 4,320 s | 7.2 GB | 9.7 to 11.5 ms wall |
| `cont`, one iteration | 30, 29, 29 s | 486 s | 7.1 GB | 1.7 ms |
| the reconstruction | 23, 23, 22 s | 374 s | 1.9 GB | 1.3 ms |
| the cone search on every row at the student's wide cones (ruling R4) | 12:05 | 16,512 s | 19.7 GB | 40.8 ms |

The in-loop cycle on trpv1 (the rows, the training, the prediction, the search of the test half,
`cont` and the reconstruction) costs 80 s of training at 180 rows, 2.5 minutes at 900, 6.8 at
4,500 and 16 at the whole half, plus 17 s of prediction and the searches, against one dense
iteration of 173 s on the same test half, and seeds nobody there.

## Skill text proposed (not edited)

A folder guide `simple-main-learn` for `src/main/learn/`:

> Use when working in SIMPLE's `src/main/learn` subsystem, the learned pose assignment: the polar
> features (`simple_polar_feats`: the student's input tensor from the polar calculator's getters,
> the conjugate wrap, the rotation by cyclic shift, the shift through the calculator's phase ramps,
> the feature cache), HHCNN (`simple_hhcnn`: the Hermitian blocks on the half ring, the
> rotation-searching template readout, the chunked trainer, the weights file) and the pose student
> (`simple_pose_student`: the model with its grid, mirror pairs, threshold, route bound and record).
> Rules: the posterior is folded over the spiral's mirror pairs before any cone is read; the student
> uses the phases and is not shift-invariant, so the augmentation stays and the shift comes from the
> search; the band is `lp`, the grid `nspace`, the model `infile` and `fname`, no other key; the
> model file is version 2, the weights file version 2; every layer must commute with conjugation
> (the suite's contract checks) and the chunked training must be bit for bit at one and sixteen
> threads; a test on simulated data needs an asymmetric phantom or the point group, `shift_sigma`
> set on a small box, and the classic generator reset before every in-process search that is
> compared bit for bit.

One line for `simple-main-strategies-parallelization`: `learn_pose_student` and `predict_poses` are
strategy routines on the caller's `(params, build)` in `simple_strategy_pose_student` with thin
private commanders; the rows come from the production fill through `get_posterior_row`. One line
for `simple-main-strategies-probabilistic`: `eul_prob_tab_neigh`'s geometric mode takes a
particle's cone from its `cone` slot about the stored direction and about `euler_mirror` of it; a
particle without the slot takes its coarse cell as before.

## The in-loop consumer, for the maintainer

The plan's consumer 2 would train a student on the rows an iteration computes anyway and seed the
neighbourhood search of the particles it did not sample. The measured cycle says: on the real sets
at 15 Å the student seeds nobody or a few percent, and its training alone costs more than the dense
iteration it would replace (trpv1: 2.5 minutes at 900 rows against 173 s); on the centred medium
assembly it seeds 46% from the whole half, which would end a sampled iteration with every particle
aligned at a third of a dense search per seeded particle, but the whole-half training's early
stopping and the real sets' broad teacher stand in the way. The proposal: not warranted as the
programs stand; a candidate only where an earlier stage centres the particles, the teacher is sharp
(a map the search already resolves) and the training is held to the cheapest design at a few
thousand rows, which costs a few minutes against an iteration's three. The sampling and
fractional-update policy is not changed by this run.

## Left open for the maintainer

- The maps of the student-only alignment (ruling R4), which on centred or large sets equal or beat
  one dense search's: whether that read is worth a route of its own.
- The training's early stopping, tuned for small sets: a tenth of the rows recovers more held-out
  seeds than the whole half.
- The side scalar, no longer read by the readout, and the shift: the student needs centred
  particles.
- The cone search's cost on wide cones (2 to 4 times the dense search); the real sets' teacher at
  15 Å (one iteration from nothing), against the plan's picture of a resolved teacher.
- A 20-seed profile of the nightly route's reports before any floor (the policy's rule); the
  acceptance values (decision 16) and the production arrangement of the halves (decision 11), both
  reserved by the rulings file.
- The in-loop consumer (above) and the multi-state extension (section 15), out of this run's scope.
- The implementation-notes row of the plan and the plan's move to `completed/` are in this phase's
  diff; the maintainer commits what is accepted.

## Commit message

```
Add the pose student: learned pose proposals for refine3D's neighbourhood search

HHCNN, a Harmonic Hermitian convolutional network over the polar Fourier
half ring, trained per dataset on the probabilistic search's own posterior
rows, proposes a direction, a cone and a confidence per particle;
refine3D refine=prob_neigh searches the cone about the direction and its
mirror. Private programs learn_pose_student and predict_poses, no new key;
projects without student records run the base path byte for byte. The
nightly route pose_student runs the split-half protocol in miniature; the
algorithm page and the completed plan record the measurements on the
Ribosembly assemblies, vcrnapol and seven real sets.
```
