# trajectory_denoise on the flex kernel: every frame a pre-image of the frames that look like it

Implementation note, 2026-10-26. **Prototype implemented and run on a first real trajectory** (written against
`7d8586b88` without a compiler at hand; the maintainer builds). First look on real data: high-resolution detail
better preserved than with `ppca`, `ppca` stronger at low resolution. Sections 1-3 are the contract and the method
as implemented, section 4 the toy behind the choices, section 5 how to run and tune it on real data, sections 6-7
validation and open points. The toy is a numpy script outside the repository (claude.ai project SIMPLE,
`claude/trajectory_denoise_flex_toy.py`).

    single_exec prg=trajectory_denoise stk=traj.mrcs smpd=0.358 pca_mode=flex nthr=16 \
        [neigs=32] [halfwin=0] [eo_pool=yes] [lp=<A>] [mskdiam=<A>] [outstk=trajectory_denoised.mrcs]

## 1. Why

`trajectory_denoise` runs `commander_ppca_denoise`; the visually best mode is `pca_mode=ppca`: Tipping-Bishop
PPCA on the real-space pixels of all frames, `neigs=200` by the commander's default, every frame replaced by
`avg + W E[z|y]`. One linear subspace has to describe every pose and state of the whole trajectory, so it cannot
hold each frame's own high-resolution detail, and a basis fitted to the noisy frames it then denoises carries
noise into every frame (with N frames < D pixels the sample eigenvalues are inflated: more components keep more
noise; at 200 components the toy's error is 44% of the raw noise energy).

`pca_mode=flex` keeps PPCA's strength, the whole trajectory informs every frame, and drops the subspace as the
image: the flex class model is only the metric that says which frames look alike; the image of a frame is a
weighted average of the data of those frames.

## 2. Contract

- In: `stk` (rows are frames in time order), `smpd`; `neigs` (largest latent rank the CV may choose, default 32),
  `halfwin` (0, the default: the whole trajectory pools; > 0: windows), `eo_pool` (default yes), optional `lp`
  (the fit band), optional `mskdiam` (the similarity mask), `nthr`. No project, no CTF (SINGLE trajectories are
  `ctf=no`).
- Out: `outstk` with the input's row count, box, sampling and row order, and `trajectory_denoise_flex.txt` (per
  frame: the Kish population of its weights and the CV's rank, tau, h).
- **No alignment.** Frames are used exactly as given: no shift, rotation, interpolation or re-centring. Every
  output row is in its input row's pose. Pose differences between frames stay in the latent (no pose nuisance
  columns), so the kernel pools only frames that look alike, pose included.
- Row t is a weighted average of input rows (only rows of t's parity with `eo_pool=yes`); row t's own weight is
  the expected weight of a frame with its clean image (section 3, step 8). A frame no other frame resembles comes
  back as itself: nothing is borrowed that the kernel does not support.
- `eo_pool=yes`: even and odd output rows share no input row. The halves of an imported stack are row parity
  (`simple_eo_partition`: a stack row without a micrograph row takes the per-index route), so the 3D half-sets of
  the denoised stack stay independent. `eo_pool=no` pools across parity (more neighbours, for viewing and 2D only).
- **No output mask.** The mask shapes the similarity only: the models, latents, noise spectrum, fit band and CV see
  soft-masked frames (`mskdiam` when given, at most the box's inscribed circle; else the inscribed circle), the
  pre-images average the frames themselves, so every output row is the whole box. Outside the mask the weights are
  the ones the particle earned, so the background is averaged over frames whose particle looks alike. Each frame's
  mean is removed before and restored after.
- Deterministic for a given input (the models are seeded per set).

## 3. Method as implemented

Files: `src/main/flex/simple_flex_traj_denoise.f90` (new, arrays in and out, outside the FLEX-PCA layering like
the cls kernel); `exec_trajectory_denoise` routes `pca_mode=flex` to `trajectory_denoise_flex` in
`single_commanders_trajectory.f90` (stack I/O); keys `halfwin`, `eo_pool` (`simple_parameters*.f90`, the UI).

1. **Planes.** Each frame minus its mean, Fourier transformed on the half-plane lattice of `flex_class_planes`
   (shells 1..nyq, `wq = 2`; `c = 1`), twice: unmasked (`yu`, what the pre-images average) and soft-masked
   (`ym`, `mask2D_soft` at radius `min(box/2 - 2, mskdiam / 2 smpd)`, what every similarity step sees).
2. **Noise.** Per shell, the median over consecutive frame pairs of `|y_t - y_t+1|^2 / 2`: model-free and robust
   to events; it includes the signal change between consecutive frames, so it errs on the high side. One spectrum
   for the trajectory; noise weights `w = 1 / s2(shell)`, the same for every frame.
3. **Fit band.** `lp` when given; otherwise up to the last shell whose spectral SNR (mean frame power over the noise
   spectrum, minus one) reaches 0.1, gaps allowed (lattice fringes sit beyond shells of weak signal), at least 6
   shells. Split into two coefficient sets A and B, Fourier blocks in a checkerboard (as the cls kernel). The mask
   correlates each coefficient with its neighbours over about box / radius pixels, so the block side grows with the
   mask: `nint(2 box / radius)`, between 4 (the inscribed circle) and 16; logged as `split block`.
4. **Sets.** `halfwin=0`: one set, the whole trajectory, every frame a candidate of every frame. `halfwin>0`: blocks
   of 10 frames, each with the frames within `halfwin` of its centre as its set (the model, CV and choice per
   block).
5. **Models.** `flex_cls_fit` on set A and on set B of the fit band, rank `min(neigs, frames - 2, coefficients/3 - 1)`,
   one nuisance column, the set mean (contrast, dose).
6. **Latents.** ML latents with the nuisance profiled out. Because c and w are the same for every frame, every
   frame's latent has the same covariance: whitened by it (unit noise per coordinate) and rotated to decreasing
   scatter. The kernel is then Euclidean:
   `w_ij = exp(-|zeta_i - zeta_j|^2_r / 2(2 + tau^2)) * exp(-(t_i - t_j)^2 / 2h^2)`, over the leading r
   directions; tau < 0: no latent term, h < 0: no temporal term.
7. **Cross-validation.** Rank r in {1, 2, 4, 8, 16, 32, 64, 128} up to the models' rank, tau in
   {0, .25, .5, 1, 2, 4, 8, 16, 32, none} (noise units), h in {2, 4, 8, 16, 32, 64, 128, none} frames, jointly.
   Set A's latents weigh the other frames (the validation frame left out, only its parity with `eo_pool=yes`) to
   predict its set-B coefficients, and vice versa; whitened error per held-out coefficient (1 is the noise floor).
   At most 1000 validation frames and 1000 held-out coefficients per set (strided). The predictions are single
   precision matrix products (`sgemm`), each thread its own block of validation frames through the whole grid, all
   temporal bandwidths of a (rank, tau) in one product. How the products run is decided once per run, with no
   environment variable: a probe compares the process CPU time of a few serial-region `sgemm` calls with their
   wall time. A serial BLAS (e.g. RHEL's `libopenblas`): blocks of validation frames, one or more per thread, each
   its own products. A threaded BLAS (an OpenMP OpenBLAS such as the bundled build, MKL, a pthreads OpenBLAS): wide
   blocks one after the other, the products from the serial region, the kernels around them in parallel; inside a
   team (`halfwin>0`) the BLAS when the build knows it stays on the calling thread there (`SIMPLE_BLAS_TEAM_SAFE`,
   see `blas_openmp_policy.md`), else the compiler's `matmul`. The log line `Trajectory flex: BLAS cpu/wall` says
   which. The optimum over the whole grid is the choice, the parity mean (tau and h
   none) included; z is its paired per-frame gain over that mean.
8. **Pre-images.** With the weights the CV scored: the fit band's set-A coefficients averaged with weights from the
   set-B model's latents and vice versa (no selection bias, so no Stein term), every other coefficient with the mean
   of the two normalised weight rows. **A frame's own weight** is not exp(0): a frame with the same clean image
   would sit at squared distance 2 per whitened direction from it (two frames' noise), so the frame gets that twin's
   expected weight, `((2 + tau^2) / (4 + tau^2))^(r/2)`; at weight 1 its own noise dominates at high rank (in the
   toy a rank-64 choice kept a third of the noise although the CV, which leaves the frame out, looked good). The
   means are matrix products over chunks of coefficients (serial BLAS: at most 128, one product per thread;
   threaded BLAS: 512, from the serial region).
9. **Readouts.** `Trajectory flex:` the band, the noise spectrum; `Trajectory flex CV by temporal bandwidth`,
   `by latent rank`, `by tau` (the curves); `Trajectory flex set=` the choice, CV, parity-mean CV, z, mean neff;
   at the end neff min/median/max. Per frame in `trajectory_denoise_flex.txt`.

Cost, whole trajectory of N frames: two model fits on half the band; per CV direction about 440 weight settings,
each a `1000 x N x 2000` single-precision product (a few minutes for N of several thousand on one node); the
pre-images three N x N weight matrices (3 x 4 N^2 bytes) and an `N x N x 2 ncoeff` product. Memory is dominated by
the frames' coefficients (8 bytes x ncoeff x N, three times: masked, unmasked, pre-images; the first two are freed
before the output is written).

## 4. Toy evidence

An fcc cluster (249 Gaussian atoms, a = 6 px, radius 15 px) in a 64 box; 600 frames in four segments tumbling at
0.1, 0.3, 0.8 and 0.05 deg/frame plus rotational diffusion, a random reorientation at each segment boundary,
halfway through the last segment a 25-atom facet cap leaves at once; tracking jitter 0.7 px rms, uncorrected; white
noise at SNR 0.2 and 1.0 (signal variance in the particle disc over noise variance). "Revisited": the last segment
retraces the first one backwards (each view perturbed by about 1 deg), no facet loss. Frames used as given; same
parity pooling; the twin weight of step 8. err = error energy over raw noise energy inside r = 24 px (raw = 1).
Global PPCA is the converged closed form of `simple_ppca` (Tipping-Bishop ML noise, posterior mean).

| err | SNR 0.2 | SNR 0.2, revisited | SNR 1.0 | SNR 1.0, revisited |
|---|---|---|---|---|
| `ppca`, 200 components | 0.444 | 0.443 | 0.455 | 0.453 |
| PPCA at its best rank (10, known only from the truth) | 0.049 | 0.047 | 0.091 | 0.081 |
| flex, window +-40 (`halfwin=40`) | 0.048 | 0.049 | 0.113 | 0.117 |
| flex, whole trajectory (`halfwin=0`) | 0.032 | 0.032 | 0.121 | 0.117 |

| SNR 1.0, standard trajectory | FRC 0.5 against the truth | FRC at the lattice fringes (raw 0.435) |
|---|---|---|
| PPCA, best rank | 5.3 px | 0.40 |
| flex, window | 3.4 px | 0.49 |
| flex, whole trajectory | 3.4 px | 0.50 |

Reading:

- Against `ppca` as it runs, 4-14 times less error. Against PPCA at its best rank: a third less at SNR 0.2; at SNR
  1.0 more total error (the global subspace is excellent at low resolution) and, on the standard trajectory,
  clearly more high-resolution detail (FRC 0.5 at 3.4 against 5.3 px); where views recur PPCA's subspace is as
  sharp (FRC 0.5 at 3.1 px against 3.1 for the window and 4.0 for the whole trajectory).
- At SNR 0.2 the whole trajectory beats the window by a third; the CV switches the temporal kernel off (h none)
  and 80% of each frame's weight comes from frames more than 40 frames away. At SNR 1.0 the window is 7% better
  where the particle tumbles fast (a local model spends its directions on the local differences) and the two tie
  when views recur. Choosing the bandwidth locally with whole-trajectory candidates did not help (it fits the
  choice to noise), and ranks up to 128 were offered: the CV took 32 or 8.
- Same-parity pooling costs 10-35% in error (half the neighbours) and removes all cross-half noise (the leak,
  the share of an output's noise variance from the other parity, is 0.42-0.50 without it).
- Without alignment the latent is what keeps the detail: at SNR 1.0 time-only smoothing halves the fringe FRC of
  the latent kernel. The ceiling is the tracking jitter: correcting the frames' shifts alone (0.25 px error)
  cut the window variant's error by a third (SNR 0.2) to almost a half (SNR 1.0); not part of this design.

Caveats: white noise, no CTF, no graphene, a near-spherical particle, 600 frames at box 64, one seed; SIMPLE's
PPCA runs 15 EM iterations from a random start where the toy uses the converged solution; the toy's latent
bandwidth is in prior units, the code's in noise units.

## 5. Running it on real data

Start with the defaults; compare against `pca_mode=ppca` on the same stack (montage of frames; then
`analysis2D_nano`). Then, guided by the log:

| look at | meaning | adjust |
|---|---|---|
| `CV by temporal bandwidth`: best at none | the whole trajectory pools | keep `halfwin=0` |
| best at small h (2-16) | the particle changes fast; time locality pays | try `halfwin` 40-100 and compare |
| `CV by latent rank`: best at the largest offered | the metric wants more directions | raise `neigs` (64) |
| chosen tau none, rank 0 | the latent carries nothing at this SNR | expected at very low SNR; check the band |
| fit band (`shells 1-k of n`) | where the metric looks | `lp=` to move it (lattice fringes in the band) |
| z small (below about 3) | barely better than the parity mean | the frames may be too noisy or too varied |
| neff per frame near 1 | unsupported frame (junk, event, lost track) | inspect those frames in `trajectory_denoise_flex.txt` |
| `eo_pool=no` | more neighbours, halves share noise | for viewing and 2D only, never for 3D |
| neighbours or carbon edges in the box drive the latent | the similarity sees more than the particle | `mskdiam=` a little above the particle (only the metric changes; the output stays unmasked) |

## 6. Validation still to do

- **Build** (done by the maintainer) and the fast gate (`scripts/check_flex_dag.py` and `scripts/check_descr.py`
  pass on the change set).
- **Unit tests** (a tester next to the module, suite in `unit_heterogeneity`): a noise-free static trajectory comes
  back unchanged; `tau` none with finite `h` gives the Gaussian temporal mean over same-parity frames (closed form);
  zero weight across parity with `eo_pool=yes`; an abrupt change with an informative latent leaves less than 1% of
  the weight across it 3 frames away; a set of fewer than 12 frames returns the frames.
- **Workflow** `test_single_workflow` (`fcc`, `wurtzite`): its trajectory is 30 orientation groups of 10 frames with
  a reorientation every 10 frames and no shifts. Run the denoise stage with `ppca` and with `flex`; add metrics
  against the clean reprojections stage 2 already writes (error over noise, FRC 0.5 against the truth); downstream
  the existing map correlation, FSC and rotation errors. `flex` becomes the default when it is better on both
  suites and downstream no worse than `ppca` within the existing margins.
- **Real data**: a GLC trajectory, `ppca` against `flex` (both `halfwin` settings), then `analysis2D_nano`. First
  run: flex keeps more high-resolution detail, `ppca` more low-resolution signal (see section 7, *low resolution*).

## 7. Open

- **The twin weight in the cls kernel route.** `flex_cls_kernel%preimage` and the representatives weigh their own
  member at exp(0) as well; for a representative's average that keeps its own noise in it. Worth the same
  correction there.
- **Tracking jitter** sets the ceiling without alignment (section 4); the tracker's sub-pixel positions as a phase
  shift would use no search and no interpolation. The owner's call.
- **Noise calibration.** The frame-difference spectrum is conservative; the CV picks the bandwidths regardless.
- **Fit band** threshold (0.1), to be settled on real data.
- **Low resolution.** `ppca`'s global subspace is strong at low resolution (as in the toy at SNR 1.0). A per-band
  choice (the CV run per resolution band, low shells with the wider bandwidths) would give each band its own
  pooling; not tried.
- **Window and whole trajectory together**: a local model as the metric with whole-trajectory candidates would take
  the best of both rows of section 4; not tried.
- **Junk axis.** Junk frames only lower their weights through their latents.
- **Cost** for long trajectories: a sparse kernel (nearest neighbours in the latent) once N reaches tens of
  thousands.
