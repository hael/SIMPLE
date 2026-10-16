# How the atomic electron tomography codes relate to SIMPLE and SINGLE

Reference note, 2026-10-14. It compares the codes of Jianwei Miao's group
(UCLA) with what SIMPLE and SINGLE do today and plan to do. The
reconstruction gets the most attention: their GENFIRE algorithm against our
PCG backend.

## Sources

**What we have read.** Peter's folder (`~/Downloads/TomoFromPeter`) holds
the full MATLAB code of the 2017 FePt study (Yang et al., Nature 542, 75):

- GENFIRE v1 and its angle refinement;
- atom tracing and species classification;
- the final position refinement;
- the data: the tilt series, the reconstruction, the traced peaks and the
  models.

On our side the sources are:

- `doc/algorithms/reconstruction.md` and
  `doc/policies/3D/reconstruct3D_pcg_policy.md`;
- `pcg_backend_overview.md`;
- `single_pcg_integration.md`;
- `species_discovery.md`;
- `atom_bfactor_model.md`;
- the code those documents point to.

**Read later.** The code of the 2021 metallic-glass study (Yang et al.,
Nature 592, 60), cloned from `AET-MetallicGlass/Supplementary-Data-Codes`
into `~/src/Supplementary-Data-Codes`. It contains:

- RESIRE, their current reconstruction algorithm;
- polynomial atom tracing;
- k-means classification with local templates;
- position refinement.

**What we have not read.** The later repository
`AET-pAET/Supplementary-Data-Codes`. Where this note describes work after
2021, it relies on the published method descriptions, not on code.

**Terms.**

- **AET (atomic electron tomography).** One nanoparticle is imaged at many
  tilts in a scanning transmission electron microscope, here with an annular
  dark-field detector (ADF-STEM). The tilt series is reconstructed to a 3D
  map in which single atoms are resolved.
- **SINGLE.** Our method for nanoparticles in a graphene liquid cell. It
  uses thousands of transmission electron microscope images of a freely
  rotating particle, each of them very noisy.
- **Gridding.** Each observed Fourier sample is spread onto the regular 3D
  Fourier grid with interpolation weights. The map is the weighted sum
  divided by the weights.
- **PCG.** Preconditioned conjugate gradients, our iterative least-squares
  reconstruction.

## 1. Two problems that end in the same place

| | AET (Miao group) | SINGLE |
| --- | --- | --- |
| Sample | one particle, in vacuum, on a support | one particle, in liquid, rotating freely |
| Images | a few tens of projections of the same particle (68 for FePt), high dose each | hundreds to a few thousand particle images, low dose each |
| Orientations | one tilt axis, about ±65°; a missing wedge of directions | random; usually no missing directions |
| Image formation | ADF-STEM: incoherent, intensity rises steeply with atomic number, no transfer function | phase contrast: contrast transfer function (CTF), weaker dependence on atomic number |
| Orientation estimate | known tilts, refined by a small grid search | unknown, estimated per image in every refinement round |
| Pre-processing | alignment and denoising of each projection (BM3D in their recent code) | motion correction, CTF estimation, 2D and 3D refinement |
| Map | 256³ at 0.3725 Å (FePt) | boxes of 120 to 200 at about 0.36 Å |

The end products are the same. Both produce a 3D map in which atoms are
resolved, an atomic model traced from it (positions, then species), and
statistics on that model: strain, coordination, chemical order.

The hard parts differ. AET has plenty of signal per image but missing
directions, so its reconstruction must fill in data that were never
measured. SINGLE has all directions but very little signal per image, and
unknown orientations. Its reconstruction must average noise correctly
across thousands of images, inside an alternation with the orientation
search.

## 2. The pipelines side by side

| Stage | AET (FePt 2017 code, and later papers) | SIMPLE / SINGLE today | Planned on our side |
| --- | --- | --- | --- |
| Reconstruction | GENFIRE (2017): Fourier gridding, then 500 iterations alternating between data, positivity and support; RESIRE (2021 on): fixed-step gradient descent on the projection misfit, with positivity (section 3) | `refine3D_nano` on the gridding backend, references filtered by ICM; PCG available in SIMPLE (`rec_backend=pcg`), not validated for SINGLE | PCG with an atomic prior inside the solve (`single_pcg_integration.md`) |
| Orientations | grid search of each Euler angle against reprojections, scored by an L1 R-factor | full pose search per image (`refine3D`, `refine=neigh`), correlation objective | unchanged |
| Atom tracing | 2017: greedy, take the brightest local maximum, fit a 3D Gaussian, check a 2 Å minimum distance, subtract the fit, repeat; 2021: every local maximum of a 3× upsampled map, each refined by a 4th-order polynomial fit, accepted brightest first at 2 Å minimum distance, nothing subtracted; manual check in both | `detect_atoms` level 1: correlation with a one-atom template, threshold, connected components, splitting, correlation validation, contact-score pruning | residual recovery of weaker atoms (`discover_species`, opt-in), which is the greedy subtraction idea applied at a second level |
| Species | intensity threshold, then global templates matched by L1 distance, manual edits, a chemistry prior near the missing wedge; later k-means and local templates | aperture intensity at one common width; Gaussian mixture with measurement noise (extreme deconvolution); BIC for the number of classes; separation gate; shared coordination slope; species list by `element=` | promotion after validation on real two-species data |
| Model refinement | positions refined against the measured projections, with one amplitude and one B per species; an interatomic potential as restraint in 2017, none in 2021 | none against the images; `_SIM.mrc` (Peng five-Gaussian kernels) becomes the next refinement reference | atomic prior in the PCG solve; per-atom B and occupancy fitted against the particle images (`atom_bfactor_model.md`) |
| Validation | multislice simulation of their own model, run through the same pipeline; R-factor against the projections; R-free in GENFIRE | synthetic element-kernel particles with ground truth; 66 real Pt maps with planted classes; half maps where they exist | image-level cross-validation (`single_pcg_integration.md`) |

## 3. Reconstruction: GENFIRE and PCG

### 3.1 What GENFIRE does

This is read from `GENFIRE_main_FePt_FinalReconstruction.m`,
`My_fill_cubicgrid_ver1.m` and `GENFIRE_iterate.m`.

1. **Gridding.**
   - Each projection is zero-padded four times (`oversamplingRatio = 4`),
     so the Fourier grid is four times finer than the image grid.
   - For every grid point within 0.05 grid units of a projection's central
     plane, the exact discrete Fourier transform of that projection is
     evaluated at the foot of the perpendicular.
   - Where several projections reach the same grid point, their values are
     averaged with inverse-distance weights.
   - Every other grid point is marked as unmeasured. So the data enter as
     nearly exact Fourier values on a sparse set of grid points. There is no
     interpolation kernel and no noise weight.
2. **Iteration.** It starts from the inverse transform of the gridded data
   and repeats 500 times:
   - set negative values to zero (positivity);
   - set values outside the support to zero;
   - Fourier transform;
   - overwrite every measured grid point with its measured value; leave the
     unmeasured ones as they are;
   - inverse transform.

   The iterate with the smallest Fourier error is kept. That error is the
   L1 distance between the measured and computed amplitudes at the measured
   points.
3. **Options used for FePt.** Measured values are enforced out to 1.5 times
   the Nyquist radius, which takes in most of the cube's corners. The
   support is a
   smoothed cube. An R-free test exists: it withholds 5% of the measured
   points in the outer shell to find the resolution the data support. The
   final FePt run does not use it.
4. **Angle refinement.** One Euler angle is refined at a time, over a grid
   of ±2.1° in 0.3° steps. Each candidate is scored by the L1 R-factor
   between a measured projection and the reprojection of a GENFIRE
   reconstruction, after a centre-of-mass shift.

### 3.2 What PCG does

This is the short version; the policy is the contract.

- **The least-squares problem.** PCG solves the weighted least-squares
  problem of the projection-slice forward model,
  `(P H P + lambda I) u = P b`:
  - `H = sum_i K_i^T K_i`, where `K_i` gathers a central section with the
    Kaiser-Bessel kernel on a 2× oversampled grid and applies the shift, the
    CTF and the noise weight `1/sigma2`;
  - `b` is the weighted backprojection;
  - `P` is a soft or hard real-space support;
  - `lambda` is a small ridge.
- **The operator.** `H` is applied as one convolution on the padded lattice
  (the Toeplitz kernel), so an iteration costs a few FFTs whatever the
  number of particles.
- **The preconditioner.** The Fourier-diagonal preconditioner is
  `1/(rho + floor)`. Here `rho` is the sampling density: the sum, over
  images, of the CTF squared over the noise power deposited at each grid
  point.
- **Iterations.** The solve starts from zero and runs two iterations in
  refinement, five for final maps.
- **Regularization.** Under ML regularization, a second, regularized map
  is the closed-form Fourier-shell shrinkage of the first, using the
  FSC-derived signal prior. SINGLE turns ML regularization off.
- **No positivity.** The solver is linear.

### 3.3 How the two relate

Written as formulas, GENFIRE is a simplified version of the problem PCG
solves, plus one constraint PCG cannot carry.

Let `M` be the indicator of the measured grid points, `m` the gridded
values, `F` the Fourier transform and `C` the projection onto
non-negative maps inside the support. One GENFIRE iteration is

    x_{n+1} = z - F^{-1} M (F z - m),        z = C(x_n).

This is one gradient step, of size one, on the gridded data misfit
`1/2 ||M (F z - m)||^2`, taken from the point `z` that satisfies the
constraints. In other words, GENFIRE is projected gradient descent.
Error reduction in phase retrieval has long been read the same way.

Now take PCG with the gridding approximation of its operator. Without
the Kaiser-Bessel coupling, `H` becomes the Fourier-diagonal density
`rho`, and `b` becomes the gridding backprojection `B`. One preconditioned
gradient step of unit size is then

    u_{n+1} = u_n - F^{-1} [rho/(rho + floor)] (F P u_n - B/rho).

Where the data are dense (`rho >> floor`), the factor in brackets is 1 and
`B/rho` is the gridded measurement. Where nothing was measured, it is 0.
That is GENFIRE's update with `M` replaced by a smooth version of itself
and without positivity. Both methods also start from the same point: the
first iterate of either is the gridding reconstruction inside the support.

PCG improves on that iteration in four ways:

- conjugate directions with optimal step lengths instead of a fixed unit
  step, so the linear part converges in far fewer iterations;
- the exact coupled operator (Kaiser-Bessel interpolation, about 3% kernel
  error) instead of the diagonal gridding approximation;
- data weights `|CTF|^2/sigma2`, so that thousands of noisy observations of
  one grid point are averaged by their reliability rather than overwritten
  by one value;
- a statistical regularizer (the FSC prior).

GENFIRE has one thing PCG lacks: **positivity**, a nonlinear constraint
that conjugate gradients cannot carry. With positivity and a tight support,
the iteration fills Fourier regions that no projection sampled, such as
the missing wedge, the gaps between few projections, and the corners of the
cube. Without positivity a support still couples measured and unmeasured
Fourier points. On the FePt tilt series PCG at first filled them harmfully,
through a preconditioner defect that is now fixed (3.5).

### 3.4 What the differences mean for SINGLE

- **Missing data against noisy data.** GENFIRE's main strength is
  recovering what was never measured. SINGLE usually measures every
  direction, many times, with very little signal each time. For SINGLE the
  noise weighting, the CTF and the averaging across images matter more, and
  GENFIRE has none of the three. Run on SINGLE data, GENFIRE would replace
  each grid point by the average of whatever samples landed within 0.05
  units, with no CTF correction.
- **Positivity as a prior for SINGLE.** The atomic potential above the
  liquid background is positive, so positivity is physically valid. At low
  signal-to-noise, though, clipping negative noise biases the map upward,
  most in the solvent and at weak atoms, which are exactly where species
  discovery looks. If we want it, it needs a projected or bound-constrained
  solver, not plain CG. On FePt it does not pay for atom detection (3.5):
  it makes the atoms rounder and the map free of negative values, but it
  does not lift more atoms above the background, and recall fell. The atomic prior of `single_pcg_integration.md`
  (penalize what is not a sum of atom kernels) is a stronger, linear
  alternative that keeps CG. It also enforces vacuum between atoms, which
  is the part of positivity that matters for atom detection.
- **Iterations and regularization.** GENFIRE runs 500 iterations and keeps
  the best one by Fourier error. That error is measured on the same data it
  enforces, so it cannot detect over-fitting; the R-free option can. PCG
  stops early on purpose: two iterations are the calibrated regime, and
  over-iteration breaks the Toeplitz approximation. Both methods rely on
  stopping as part of their regularization.
- **Cost.** GENFIRE's exact gridding costs about the number of projections,
  times the grid points near each plane, times the pixels per projection.
  That is affordable for 68 projections and impossible for thousands of
  images. PCG reads every image once, and its iterations do not depend on
  the number of images.
- **RESIRE.** Their current algorithm (Pham et al., Sci. Rep. 2023),
  read from `2_RESIRE_package/src/@RESIRE_Reconstructor/reconstruct.m`.
  It minimizes the squared difference between the measured and computed
  projections by gradient descent with a fixed step:
  - **Computing projections.** The volume is padded four times and Fourier
    transformed. The central plane of each projection is interpolated with
    splines, then transformed back to real space.
  - **The gradient.** It is the real-space back-projection of the computed
    projections, minus the back-projection of the measured ones (computed
    once, before the loop).
  - **The step.** The step is fixed, scaled by the number of projections
    times the box size. Negative values are set to zero after every step.
    The default run is 200 iterations.
  - **What is not there.** The loop has no support constraint (a support
    can only serve as the starting model) and no noise weights. The forward
    and back projections use different interpolations, so they are not
    exact adjoints.

  This is PCG's objective without the weights, solved by fixed-step
  gradient descent with positivity instead of conjugate gradients. It is
  also GENFIRE with the misfit measured on the projections rather than on
  gridded Fourier points. Their field has moved toward our formulation,
  but it keeps positivity and has no noise model.

### 3.5 A direct comparison we can run

SIMPLE can reconstruct a tilt series. The projections become particles
with known orientations (after converting their Euler convention),
`ctf=no` and the correlation objective. The FePt
data allow a like-for-like comparison: `reconstruct3D rec_backend=pcg` on
their projections against their GENFIRE volume, which is in the folder.

It needs the denoised projections and refined angles from the group's FePt
page; the folder holds only the raw 300² projections. The comparison would
show three things:

- how much of their map's quality comes from positivity in the missing
  wedge;
- what PCG's support does there without positivity;
- whether our atom tracing and species calls change between the two maps.

**First run, 2026-10-15.** It used the denoised projections with the
initial single-axis angles (`for_simple/run_fept_recon_only.sh`), and the
frame and hand came out right.

- **Gridding works.** The half-map FSC reaches 0.5 at 3.4 Å.
- **PCG did not.** It reached 0.5 at 6.6 Å after one iteration and 9.2 Å
  after two, with a checkerboard spike at the map centre, 13% of its power
  beyond the Nyquist sphere and a halo in the missing wedge. One defect
  caused all three: the right-hand side is cut by the solve support before
  the Fourier divide, and in voxels no projection measured the
  preconditioner was only the floor, about a hundredfold gain. The cause,
  the fix and its test are in
  `doc/refactoring_notes/completed/pcg_sparse_sampling_fix.md`.

**After the fix, same day.** PCG reaches 0.5 at 3.41 Å after two iterations
and 3.31 Å after five, and correlates better than gridding with a numpy
filtered back-projection (0.919 and 0.911 against 0.896). Level-1 atom
detection, against their final model, finds:

| Map | Recall | Precision |
| --- | --- | --- |
| SIMPLE gridding | 0.25 | 0.99 |
| SIMPLE PCG, 2 iterations | 0.26 | 0.99 |
| SIMPLE PCG, 5 iterations | 0.35 | 0.99 |
| Their GENFIRE map | 0.41 | 1.00 |

With species discovery on the half maps, PCG at five iterations finds 54%
of the atoms at 98% precision (gridding 41%). The species calls are wrong:
the detected atoms are almost all Pt, since few Fe atoms rise above the
background in any of the maps, and discovery splits the Pt into two groups.

**Positivity.** Their atoms look rounder, which raised the question whether
positivity is the reason. A numpy GENFIRE-style iteration on our data and
nominal angles (measured Fourier points put back each iteration, a sphere
as support), with and without clipping negative values:

- **It rounds the atoms.** The mean Pt atom is 2.02 times longer along the
  missing-wedge direction than across without positivity, 1.90 after 200
  iterations with it. Their map: 1.89.
- **It does not separate atoms from the background.** The share of Pt
  atoms brighter than the 99th percentile of the background between atoms
  stays at 85% (theirs 100%).
- **It lowers detection.** Level-1 recall is 0.30 without positivity and
  0.24 with it.

So positivity explains the look of their map, not its better detection.

**Refined angles.** `refine3D refine=cont` refines each projection's three
rotation components and two shifts against the PCG map of the nominal tilts
(`for_simple/run_fept_refine_cont.sh`). It converged in two iterations and
moved the projections by 0.31 degrees on average (largest 1.1 degrees) and
0.2 pixels, with no rotation of the whole frame. PCG at five iterations from
the refined poses raises level-1 recall from 0.35 to 0.39, at precision
0.99 (their map: 0.41), and recall with discovery from 0.54 to 0.62. The
share of Pt atoms above the background rises from 66% to 71% (theirs 90%,
smoothed measure). Angle refinement in our own pipeline therefore closes
most of the detection gap to their map. The SINGLE workflow `refine_tomo` packages
these steps, with a density mask as the PCG support and, by default, two
iterations at 2 A: longer refinement stretches the tilt range and rotates the
frame, which the data barely constrain; `import_tilt_series` reads their angle
convention and `scripts/aet_mat2simple.py` their MATLAB files
(`doc/refactoring_notes/planned/refine_tomo.md`).

## 4. Atom tracing

| | AET | `detect_atoms` |
| --- | --- | --- |
| Candidate atoms | every local maximum above a threshold, brightest first | voxels of the template-correlation map above a threshold, chosen by how well the resulting atoms simulate the map, grouped into connected components |
| Separating close atoms | each local maximum is its own atom; a 2 Å minimum distance (in 2017 the fitted Gaussian of each accepted atom was also subtracted before the next maximum was taken) | merged components are split, validated against the correlation map |
| Position | 2017: centre of an anisotropic 3D Gaussian fitted in a 5 × 5 × 5 box; 2021: maximum of a 4th-order polynomial fitted on the 3× upsampled map | centre of mass of the component, after the splitting |
| False atoms | a "non-atom" class in the classification, then manual removal | correlation validation and contact-score pruning (weak, poorly coordinated atoms are removed on purpose) |
| Weak atoms | found by the same greedy loop as strong ones | level 1 finds the strong class; `discover_species` recovers weaker atoms from the residual at calibrated noise thresholds, with a neighbour gate |

Their tracing needs no segmentation and treats every atom alike. Ours
segments first. The residual recovery of `discover_species` is their 2017
idea applied only below the strong class, with thresholds set against a
measured noise level instead of a fixed intensity threshold and a manual
check. In the papers we looked at, hand edits are part of their pipeline.
Ours has none, and that is a design constraint we keep.

The FePt map shows what segmenting first costs on a dense two-species
particle (2026-10-15). With `element=Pt`, level 1 found 57% of the Pt atoms
and none of the Fe. Every atom it found is real, within 0.14 Å of their
refined position.

The atoms it missed are not merged into components with found ones. They
lie below the threshold. The search chose that threshold because lower ones
merge neighbouring atoms: the correlation between the map and the
simulated atoms falls from 0.50 to 0.28 as the threshold drops. A tracer
built on local maxima does not tie the threshold to merging.

## 5. Species

Their methods across papers:

- **FePt (2017).** A threshold on the 3 × 3 × 3-voxel intensity; global
  5 × 5 × 5 templates matched by L1 distance and iterated; a test against
  zero for non-atoms; manual edits; reclassification of Fe atoms near the
  missing wedge with a separate threshold.
- **The high-entropy metallic glass (2021) and nucleation (2025) studies.**
  K-means on integrated intensity into three types and non-atoms, then
  local templates: each atom goes to the type whose mean within 10 Å it is
  closest to.
- **The PtNi catalysts (KAIST group, 2025).** K-means into non-atom, Ni and
  Pt, with an Otsu boundary.

None of these models how an atom's environment changes its intensity. On
the FePt volume (`species_discovery_status.md`):

- intensity rises by 5.2% per neighbour for both species;
- the classes are separated by `D = 2.6`;
- the share of Fe at the surface is 65 to 68% in their model, 86% from a
  plain mixture fit and 20% from the same fit with our coordination term.

Local templates do not remove the trend (87% at the surface), because the
10 Å neighbourhood of a surface atom is mostly bulk.

Ours: a mixture with an explicit measurement noise, a criterion for the
number of classes, a separation gate, and one shared coordination slope.
The comparison says the environment model is the decisive choice for
surface species. Neither group can settle it from map intensity alone.

## 6. Refining the model against the data

Their final refinement (`XYZ_Refinement_Main.m`):

- **The computed data.** Each atom contributes, in Fourier space, a common
  form factor (the average of the Kirkland Fe and Pt factors) times a
  per-species amplitude and a per-species Debye-Waller term. The
  amplitudes are 0.56 and 1.18; the B factors are 21.9 and 20.1 Å² in the
  usual convention. The sum is compared with the transforms of the 68
  measured projections, inside a resolution sphere, after a least-squares
  scale.
- **The steps.** Positions move along the gradient of the squared
  difference plus a weighted gradient of an embedded-atom potential for
  FePt. Each step has a fixed total length. After 150 steps, the
  configuration with the smallest data error is kept.
- **What stays fixed.** Species and B are not refined per atom, and there
  is no occupancy.

On our side:

- **Today.** Positions come from the map, and `_SIM.mrc` (Peng's
  five-Gaussian electron scattering factors per element, with a B) only
  seeds the next refinement round.
- **The atomic prior** (`single_pcg_integration.md`). It puts the model
  inside the reconstruction: `lambda_A (I - P_G)` penalizes density
  outside the span of the atom kernels, with amplitudes free and scale
  never compared. Its strong-prior limit is their refinement's model (a
  sum of atoms fitted to the data), reached through the image data rather
  than through a separate refinement.
- **The B-factor note's preferred route.** It fits a smooth per-atom B and
  occupancy profile against the particle images. That is their
  refinement's objective with the two per-atom quantities they lack.

Their B factors are of the same size as the total B we measured in the
core of real SINGLE Pt maps (about 21 Å²).

A related line of work, local-orbital tomography (Rong Yu's group,
Tsinghua), makes atoms the unknowns of the reconstruction, with a
depth-dependent probe. That is the far end of the atomic-prior ladder.

## 7. Validation practice

- **Their checks.** They simulate their own final model with a multislice
  code, run the simulated projections through the same pipeline, and report
  the fraction of atoms recovered with the right position and species. They
  also report the R-factor of the model against the measured projections.
  These checks test self-consistency, not truth.
- **The independent critique.** Busch et al. (Nature, 2026) argue, from
  simulations, that species identification becomes very uncertain once
  atom-peak intensities overlap the background. The group's reply in the
  same issue reports better recovery, but also on simulated data, so the
  question stays open for real data.
- **Ours.** Ground truth exists only for synthetic particles (element
  kernels, realistic noise). On real maps we test single-species maps for
  false species and plant weaker classes. The FePt volume is now our first
  real two-species map, with their model as a reference, not as truth.

## 8. What we take from them

- **The FePt data,** as the first real two-species test of
  `detect_atoms`, with and without the coordination term
  (`species_discovery_status.md`). It is also a reconstruction test: PCG
  against GENFIRE on the same projections (3.5), once the denoised
  projections and angles are downloaded.
- **Positivity, considered as a prior for SINGLE.** It has a known
  low-signal bias, and the planned atomic prior covers its useful part
  within a linear solve.
- **Refinement against the images, not the map,** for positions and,
  unlike theirs, per-atom B and occupancy. It is the one place where
  occupancy and species might be told apart.
- **An R-free style hold-out for every iterative or fitted step.** It is
  their GENFIRE option and our planned image-level cross-validation, in
  place of a fit criterion measured on the data being fitted.
- **What we keep:** no manual edits, explicit noise models, and an
  environment term reported with every surface species call.
