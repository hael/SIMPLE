# The Pose Student

## Problem

Given a set of particle images whose poses have been estimated by the
probabilistic 3D search ([refine3D](refine3d.md), `refine=prob`), predict for a
particle the search has never seen a projection *direction*, a *cone* about it
that holds most of the posterior, and a *confidence*, so that the search can
align that particle inside the cone instead of over the whole sphere. The
student proposes; the search decides. It never writes an in-plane rotation or
a shift: those come from the neighbourhood search it hands the cone to, with
the search's own machinery (the polar shift map, the refinement, the local
refinement).

The quantity estimated is a distribution over the *mirror pairs* of the
student's grid, modulo the point group, for one state. The grid is the
search's reference spiral of `nspace` (1,000) directions, whose first half
holds the unmirrored directions of the asymmetric unit and whose second half
their mirrors (`build_refspiral`; the pairs and their frame angles are what
`oris%get_mirror_pairs` returns, the pairing the opt-in `mirr_proj=yes`
reference extraction uses). A direction and its mirror view are related by the
Euler angles (`euler_mirror`), and the projection of one is the reflected
image of the other: the student classifies the 500 pairs and leaves the
choice within a pair to the search, which resolves it with a cross-correlation.
The pairing is the three-dimensional Friedel symmetry seen in projection: a
real volume's transform obeys F(-k) = conj F(k), so the central section for
the antipodal direction is the same plane in a reflected in-plane frame. A
reflection of an image is itself a reflection of its Fourier plane, θ → -θ,
not a conjugation; conjugation enters only because the half ring stores
[0, π), and Friedel brings -θ back as F(k, -θ) = conj F(k, π - θ), so on the
half ring a reflection is an index reversal with a conjugation, and in the
angular spectrum of the full ring the reflection θ → π - θ maps F(n) to
conj F(n). For the cyclic point groups, c1 included, the unmirrored half is the
northern hemisphere (polar angle below 90°); for the dihedral and platonic
groups it is cut in the azimuth instead.
Everything continuous in the pose (the in-plane angle, the shift, the
direction within the cone) is left to the search as well.

## Model

**The input** is the particle's polar Fourier half ring exactly as the search
holds it: `pftsz` angles over half a turn, the rings of the band `kfromto(1)`
to `box smpd / lp` (about 25 rings at a 15 Å band for a box of 128 pixels at
3 Å), the other half of the ring being the Friedel mate. The ring is multiplied
by the particle's CTF (the cross term of a CTF-free model with the observation,
as the search's objective and flex's Gram both correlate), each shell scaled by
the square root of the shell weight `k / sigma2(k)` the search scores with, and
the band scaled to unit weighted power; the ring mean of the CTF squared
travels beside the tensor. Nothing is restored (no phase flip, no Wiener
filter) and the phases are kept. A shift is applied, when the training needs
it, through the search's own phase ramps; an in-plane rotation by one grid
step is a cyclic shift of the angular axis whose wrap conjugates what passes
the ring's end.

**The network's half ring** is the calculator's capped at the band. A ring of
Fourier radius k holds angular frequencies up to about k, so a half ring needs
π k_max angles. The calculator's `pftsz` follows the mask radius instead and
oversamples a low-pass band: 220 angles on bgal (256 pixels at 1.275 Å) where
its 15 Å band of 22 rings needs 69. The network takes at most
`NET_ANGLE_OVERSAMPLING` (2) times π k_max angles, even
(`polar_feats_net_pftsz`). Each ring is resampled exactly onto them: the
Friedel-extended ring's angular harmonics |j| < `pftsz_net` are kept, as two
precomputed real matrices on the real and imaginary channels. The scale keeps
the band's unit weighted norm. The teacher, the calculator and the model's
`pftsz` stay the calculator's, so `predict_poses` checks the geometry as
before. Where the cap is not below `pftsz`, as on the medium assembly (128
pixels at 3 Å, 160 angles, cap 164), nothing changes.

Measured on 2026-11-11 against six base seeds, three seeds per cap, at 24
threads on the Dell:
- bgal at twice π k_max (140 angles): within 6° 0.876 against the base's
  0.869 (range 0.853 to 0.880), truth in the cone 0.959 against 0.957, seeded
  0.426 against 0.417. The designs train in 15.8 s against 32.8 s, and
  prediction's forward pass takes 1.33 s against 2.37 s.
- 1.5 times (104 angles): the reads hold on bgal, but T1000's seeded fraction
  falls to the edge of its range (0.216 against 0.215 to 0.272).
- 1 time (bgal 70, T1000 82): the student loses on both sets (T1000 within
  6° 0.756 against 0.784, seeded 0.191 against 0.242). This is the loss the
  low-pass crop of `box_crop` showed (`doc/implementation_notes/rejected/pose_student_working_grid.md`),
  so it came from the angles.

**The network is HHCNN, the Harmonic Hermitian convolutional neural
network**: Hermitian because every layer commutes with conjugation on the half
ring with the conjugate wrap, so the network is a model of the Fourier
transform of a real image; harmonic because the readout's class templates live
in the angular harmonic domain of the polar transform and the rotation search
is the cross-correlation over the cyclic shifts through the angular Fourier
transform; convolutional on the polar Fourier representation. It is not
Worrall's Harmonic Networks (circular-harmonic steerable filters): the
convolution kernels here are ordinary real kernels in the ring-and-angle
domain. Three blocks of a 3 by 3 circular convolution along the angle and the ring with real
kernels applied alike to the real and imaginary channels (the only linear map
of a complex signal that commutes with conjugation), a bias on the modulus and
the gate `z relu(|z| + b) / |z|` that mixes the two parts only through the
conjugation-invariant modulus, and an average over ring pairs. The kernels
are symmetric along the angle (their two angle offsets share one weight), so a
reflected input gives reflected feature maps; training keeps them symmetric
bit for bit, the gradients of the two offsets being summed onto the shared
weight. The first
block's filters are CTF-free and each response is divided by the filter's
energy under the particle's ring-averaged CTF squared and the shell weights,
so that the CTF's envelope never reaches the deeper blocks as structure and a
shell near a zero enters at no weight.

**The readout is a learned matched-filter bank with a rotation search.** A
pointwise map with real weights reduces the last block's channels to eight, an
average over ring bands reduces the rings to eight, and the band-limited
angular spectrum of the full ring is taken to order 16: its even orders from
the real part of the half ring and its odd orders from the imaginary part,
which is what the Hermitian extension gives. Every direction class owns a
template in that domain (eight by eight by sixteen complex values, the
order-0 term real), initialised as the spectrum of the reprojection of the
pair's unmirrored direction through the untrained network and then trained
with the rest. The score of a class is the maximum, over the ring's rotation
grid and over the map and its reflection, of the real Hermitian inner product
of the particle's spectrum with the template,

```text
S_c(r) = (1/N) [ Re P_c(0) + 2 sum_{m=1..16} Re( P_c(m) e^{2 pi i m r / N} ) ],
P_c(m) = sum_{channels, bands} F(m) conj(T_c(m)),
```

with `N` the rotation grid's length, and for the reflected map `F(m)` replaced
by its conjugate (the reflection of the ring in the angular-harmonic domain);
the logits are a learned scale times the maxima, the posterior is their
softmax, and the gradient flows through the maximising rotation and, where the
reflection won, through the conjugate spectrum. Rotation equivariance is exact
on the search's rotation grid: rotating the particle by one step leaves every
class score unchanged and moves the maximising rotation by one step (backwards
for a class the reflection scored). With the symmetric kernels the network is
invariant to the mirror on that grid: a particle and its reflection, mirrored
by the calculator's `mirror_ref_bank` at any frame angle on the rotation grid,
score alike against every class. The readout keeps the angular structure that
tells views apart; an earlier head that pooled rotation-invariant statistics
(the angular mean of the real part, the mean modulus, harmonic magnitudes)
discarded it, which is why it was replaced.

## Training

**The teacher** is the probabilistic search itself: for every training
particle, the row of the dense probability table at the student's grid and
band, filled as the production search fills it against the training set's map,
each direction scored at its drawn in-plane rotation with the best candidates
refined in shift. The soft target is the row's tempered posterior
`t_c = exp[-tau (L_c - L_min)] / sum_j exp[-tau (L_j - L_min)]`, `tau` the
search's tempering factor times the particle's own temperature; the whole row
is kept, so a broad posterior stays broad, and it is folded onto the mirror
pairs (a pair's target is the sum of its two directions') before the loss. The
teacher scores all 1,000 directions; only the student's classes are pairs. The
shell weights the rows were whitened with are frozen in the model.

**The loss** is the cross-entropy against those soft rows, minimised with Adam
(rate 1e-3, decoupled weight decay 1e-5) on mini-batches of at most 256, every
row shifted each epoch by a keyed draw in a disk of three times the search's
centring width (the training's shift tolerance), for at most 30 epochs and
stopped after 3 epochs without a held-out improvement (the best epoch comes at
3 to 7 on the fixed sets, its improvements at most 2 epochs apart); the best
epoch's weights are kept.

**The training step** is built so that the same weights come out at any
thread count, bit for bit. The mini-batch is cut into chunks set by its size
alone (about 96 chunks), computed in parallel and reduced in chunk order over
slices of the weights, also in parallel, with the Adam update in the same pass.
The blocks' convolutions are im2col and a matrix product; blocks 2 and 3 build
their im2col matrix in tiles of 128 columns, which stay in a core's cache. The
readout's class products, its rotation scores and their gradients are matrix
products too; the rotation search scores every rotation of the grid (a
coarse-to-fine search was measured and changed the winning rotation too often
for its gain), and the templates' gradient is one product over the whole batch.
Every matrix product takes the same arithmetic whatever the team: `sgemm`
where the build's BLAS stays on the calling thread inside OpenMP teams
(`SIMPLE_BLAS_TEAM_SAFE`). Otherwise the build and the shape choose between
explicit loops down the columns of the result and the `matmul` intrinsic,
gfortran's single-threaded product, which is compiled into its runtime and so
runs optimised even in a Debug build. On the Dell's Arrow Lake core, one
thread, the optimised loops ran at 30 to 54 GFLOP/s on the network's products
against `matmul`'s 6 to 43. The exception is a transposed product of fewer than
16 rows (the readout's four): there the loops ran at 6 and `matmul` at 37. In a
Debug build the loops ran at 4 to 5 and `matmul` at 9 to 63. So a Release build
takes the loops except on those few-row products, and a Debug build takes
`matmul` everywhere, transposing an operand first where `matmul` of a
transposed section would run at a sixth of the speed.

**Held-out micrographs** (every tenth of the selected particles' stacks in
order of their numbers when they lie on ten or more, otherwise every tenth
selected particle; fewer than ten selected particles are refused before
training, so that neither cohort is empty) serve three decisions, and the test
half is never touched before the model is frozen: early stopping on the
held-out loss; the choice of the design, where the cheapest (16 channels) is
kept unless 32 or then 64 channels recover clearly more held-out seeds (the
paired difference over micrographs, resampled a thousand times by micrograph;
kept when the 5% quantile of the gain is above zero; 64 channels are tried
only from 10,000 training rows: below, they rarely beat 32 and cost two and a
half times as much); and the confidence threshold.

**A seed is recovered** when the teacher's direction, or its mirror partner,
lies inside the student's cone, the cone being within the cap of 64 directions
and within the route's bound (below). The threshold is the confidence of the
longest confidence-sorted prefix of held-out particles in which at least 95%
of the seeds are recovered, a prefix ending only where the confidence drops
(prediction takes every particle at or above the threshold, so a prefix inside
a group of equal confidences would seed the whole group); above every
confidence when no prefix qualifies, so that a student that cannot be trusted
seeds nobody. The threshold and the bound are written so that they read back
bit for bit, and learn_pose_student reads the model back and checks that it
seeds the calibrated held-out particles.

## Prediction and the hand-off

**The readout's pose.** The readout already estimates the in-plane rotation
and the mirror: it scores each class at every rotation of the polar grid,
against the template and its reflection, and only the best score enters the
softmax. Read out for the mode's pair, the best rotation `ang` gives a full
pose with no search: the unmirrored direction with `e3 = 360 - ang` (the
search's convention for its rotation index), or, when the reflection scored,
the mirror partner with `e3 = ang + psi`, `psi` the pair's frame angle (the
reflection `mirror_ref_bank` applies; the HHCNN tester pins it at 0, 75 and
180 degrees). The template's frame is the unmirrored direction's reprojection
at `e3 = 0`, where its training starts. The network does not estimate the
shift: it is trained to tolerate it. Given the map (`vol1`), `predict_poses`
registers each particle's shift against the map's reprojection at the
read-out pose, its direction and rotation fixed: a 5 x 5 grid over 25 Å
about the stored shift, then L-BFGS-B, the search's own fixed-rotation shift
search (`pftc_shsrch_grad`), with the references and sigma2 prepared as
learn_pose_student prepares them; shared memory only. The table's last five
columns are the pose and the shift (the stored shift without `vol1`); the
route and the records do not use them.

The posterior is **over mirror pairs**: each pair's mass is placed on its
unmirrored direction of the grid, the mode and the smallest set holding 99% of
the mass (the cone) are taken over the pairs, and the cone's radius is that
set's spread about the mode, a member's distance being the smaller of its own
and its mirror's. The direction written is the mode pair's unmirrored
direction; which of the pair the particle shows is the search's to decide.

**The route** sends a particle to the neighbourhood search when its cone is
within the cap, its radius is within the bound (the dense search's cell:
`prob_athres` or the grid's spacing, whichever is larger, recorded in the model
with the threshold calibrated against it) and its confidence is at or above
the threshold; every other particle takes the dense search. The seeded
particle's record carries the direction as its orientation and, in numbered
slots of the particle's record, the cone's radius, the confidence and the
route. The program writes a part file for a distributed caller to merge, or
its own project copy, never the caller's project. A part is given its range as
the distributed masters hand it to every part (`part`, `nparts`, `fromp`,
`top`, with `outfile`); its tables carry its number, and a range without a
selected particle writes its records unchanged. Prediction refuses a model of
another geometry than the project's: the half ring, the box, the crop
(`box_crop`) and the CTF treatment exactly, the sampling (`smpd`) within
1e-4 Å and the mask (`mskdiam`) within 0.01 Å, the header's decimals; it does
not resample.

**The search** (`refine3D refine=prob_neigh prob_neigh_mode=geom`) scores, for
a particle with a cone, every reference within the cone about the stored
direction and every reference within the same cone about its mirror (`euler_mirror`), finds
the in-plane rotation and the shift as it does for any neighbourhood, and
takes the better pose as it takes any; a particle without a cone runs the base
`geom` search, the directions within `athres_cont` of its stored direction
(the coarse cell until 2026-11-09), so a project without student records runs
the base code path byte for byte.

## What the reads found

The split-half protocol (view-balanced halves, the training half's map as the
references, the student trained on one half and judged on the other against
the standard machinery's alignment and, on simulated sets, the truth) gave
these results on the Dell in November 2026; the measurements are in the
implementation note's Progress record.

- On centred particles of the simulated Ribosembly medium assembly at 15 Å the
  student seeds 46% of the test half with cones about 8° wide, those particles
  end where the dense search puts them (99% within 6° of the truth), the
  prediction costs 0.56 ms per particle against 2.2 ms for the coarse fill it
  is trained on, and the cone search costs 0.6 of the dense search. On the
  easy simulated set (vcrnapol, SNR 0.1) it seeds 86%.
- The student is not shift-invariant: a 2 px shift changes its posterior by a
  total variation near 0.7, and on uncentred particles (no prepared shift) it
  seeds about 1%. The shift is the open question of the design.
- On the real sets at 15 Å (trpv1, beta-galactosidase, streptavidin, the
  flagellar export gate FlhB, two CESPED entries) one dense search from nothing
  is itself a broad teacher (6 to 69% of its directions within 6° of the
  reference alignment), the student's confident cones are few (0 to 4% seeded),
  but reconstructing from its directions alone, with the cone search for the
  in-plane rotation and shift and one local-refinement iteration, gives maps as
  good as or better than that dense search's on centred or large sets.
- A tenth of the training rows can recover more held-out seeds than the whole
  half: the training budget's early stopping is tuned for small sets.

- Classifying the 500 mirror pairs instead of the 1,000 directions, with the
  mirror left to the search, costs nothing and gains with more training data.
  On the centred test half of the medium assembly at 15 Å, the student trained
  on the whole training half puts 94.0% of the particles within 6° of the truth
  against 91.5% for the 1,000-class student, with 99.5% within 12°; its cones
  are tighter (a 90th percentile of 14 directions against 42), and 97% of the
  particles clear the confidence of 0.99 with the truth inside the cone for
  99.7% of them. From 1,000 training particles it puts 85.5% within 6° against
  81.7%; from 200 the two are level (about 70%). One seed per point.

## Rationale

The student is a classifier over a discrete grid because the posterior over
directions is multimodal at low signal-to-noise and the references are
discrete; a softmax over the grid holds competing directions with their
masses, which a regressed rotation cannot. Its classes are mirror pairs
because the mirror is a relation of the Euler angles that the search resolves
for free, with a cross-correlation of the two candidates: a network that had to
tell a view from its mirror would spend its capacity on a distinction the
search makes exactly, and would learn half as much data per class. It hands over a cone rather than a
pose because the search already has the better continuous estimator for the
in-plane angle, the shift and the direction within the cone, with the
references, the shift map and the refinement. It is trained on the search's own
posterior rows rather than on hard labels so that it inherits the search's
uncertainty and never claims more than the teacher resolves: where the
teacher is broad, the student seeds nobody and the search runs as it would
have.

## Implementation

`src/main/learn/`: `simple_polar_feats` (the input tensor, the augmentation,
a disk cache of the rows' tensors), `simple_hhcnn` (HHCNN: the Hermitian blocks,
the template readout, the trainer, the weights file), `simple_pose_student`
(the model with its grid, its mirror pairs, threshold and record; a text
header with the weights beside it, model version 4: the header names the
weights by their file name alone, resolved against the header's directory and
never the working directory, and the reader refuses weights whose shape, seed
or weight count differ from the header's; the header's records may be of any
length, as the shell weights of a wide band need). The programs
`learn_pose_student` and `predict_poses` are
strategy routines on the caller's `(params, build)` with thin private
commanders; the band is `lp`, the grid `nspace`, the model `infile` and
`fname`; no other key. The per-particle cone of the neighbourhood search lives
in `simple_eul_prob_tab_neigh`. The design, the rulings that shaped it and
every measurement are in
`doc/implementation_notes/completed/learned_pose_state_assignment.md`.

## Validation and benchmarks

**`validate_pose_student`** (simple_exec, developer visibility) takes solve3D's
inputs on a project with a 2D solution; its `nsample` is the training set: a
class-balanced greedy sample of the selected classes, best 2D scores first.
It runs solve3D on that sample with `nsample` at the sample's size, so every
iteration searches every training particle (no fractional update, no
trailing reconstruction). It trains the student on solve3D's poses and map at
15 Å, with a shift disk of three times 2 Å since every particle is centred by
the 2D or 3D solution, and predicts every selected particle: the training
rows at the teacher's poses, the others at their 2D shifts. The report
(`validate_pose_student_report.txt`, also in the log) gives:

- the model's threshold and held-out recovery;
- the seeded and capped fractions, the confidence and the cone radius of all,
  training and unseen particles;
- the student's mode against the teacher on the training particles;
- class consistency: the distance of every particle's direction to its 2D
  class's teacher medoid, the teacher's own spread being the baseline;
- the training rows' shift change from the 2D solution to the teacher's,
  against the training's shift disk;
- the readout's full pose against the teacher's on the training particles
  (the direction with the readout's mirror, and the in-plane angle);
- the registered shifts against the teacher's, beside the 2D shifts;
- six maps by reconstruct3D from poses with no orientation search: the unseen
  particles at the student's poses with their 2D shifts and with the
  registered shifts (their FSC against the teacher's map is cross-validated:
  the two share no particle), the training particles likewise, the training
  particles at the student's poses and the teacher's shifts, and the
  teacher's own poses as the control; each with its FSC against the teacher's
  map and between its halves;
- the phases' wall times with the per-particle costs.

The phases' own records are `BENCH/VALIDATE_POSE_STUDENT_BENCH.txt`, solve3D's
`BENCH/` in its directory, `learn/BENCH/POSE_STUDENT_LEARN_BENCH.txt` and
`predict/BENCH/POSE_STUDENT_PREDICT_BENCH.txt`.

**`POSE_STUDENT_LEARN`** splits learn_pose_student into the rows (setup,
particles, the teacher's table fill with its fill profile, the features) and,
per design, the training into batch assembly, gradients (the parallel
region), reduction with the Adam steps and the held-out evaluation, with the
steps, the samples per second, the held-out prediction and the design choice.
**`POSE_STUDENT_PREDICT`** splits predict_poses into its phases and the
per-particle steps in thread-seconds (features, forward pass, cone, the
reference error, ruling R2's shift check) with the team's utilisation and load
balance. **`test=pose_student_bench`** (simple_test_exec, no CTest entry) times
the HHCNN kernels on synthetic tensors per ring count and width: the forward
pass split into the blocks and the readout, the backward pass, a training step
of 256 samples at one thread and at the team, the team's prediction throughput
and the step's identity across threads; the arithmetic rates of the blocks'
and the readout's products, the chunks per step and the threads they keep
busy; and the blocks' other kernels (im2col, the gate, the pooling, their
backward passes) at one thread (`BENCH/POSE_STUDENT_KERNELS_BENCH.txt`).
