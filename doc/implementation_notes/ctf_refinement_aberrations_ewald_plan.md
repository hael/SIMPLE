# Defocus refinement, higher-order aberrations and the Ewald sphere: survey and plan

## Status

**Planning only (2026-10-01). No code is proposed for immediate implementation.**

The order of work in section 6 was agreed by Hans on 2026-10-01. The first
item is per-particle defocus and astigmatism refinement with analytic
derivatives, driven by SIMPLE's own L-BFGS-B optimizer and built on the same
Cartesian framework as continuous pose polishing. That framework is under
heavy refactoring and testing
(`doc/refactoring_notes/planned/pose_cont_refactoring.md`), so nothing here
starts before it has landed. This note records the findings and the plan; it
is not a design.

Sources read: RELION 5.1.0 (`e18cc771`), cisTEM master (`455de9a1`), specter
main (`6832b98`), SIMPLE (`c961fb432`). cryoSPARC is closed source and is
described from its public guide and forum. specter was read, not run.

## 1. Why this matters now

`refine3D_auto` reproduces published resolutions from the same particles on
noise-limited targets and then stops improving. What is left is the data
model: each particle carries a CTF taken from a smooth per-micrograph defocus
surface, there is no beam tilt or higher-order aberration model, the
reconstruction assumes a flat Ewald sphere, and particles are not polished.
Polishing has been started separately. This note covers the CTF side.

## 2. What the other packages do

| | RELION (`relion_ctf_refine`) | cisTEM (`refine_ctf`, `refine3d`) | cryoSPARC | specter (Ghostbuster) |
| --- | --- | --- | --- | --- |
| Per-particle defocus | L-BFGS with analytic gradient on `sum w abs(Y - CTF P)^2`, all particles of a micrograph jointly; phase shift, defocus, astigmatism, Cs and B/scale each per particle, per micrograph or fixed (`--fit_mode`, default `fpmfm`: defocus per particle, astigmatism per micrograph); optional brute-force 1D pre-scan over +-2000 A | Brute-force 1D scan, +-500 A in 50 A steps, astigmatism fixed, scored by the Frealign weighted correlation at fixed pose | Per-particle defocus search (Local CTF Refinement), also on the fly inside refinement | One global scalar offset, declared trainable but not connected to the forward graph |
| Odd aberrations | Per optics group: `sum CTF conj(P) Y` and `sum CTF^2 abs(P)^2` accumulated per Fourier pixel over all micrographs, then beam tilt or odd Zernikes fitted (`--odd_aberr_max_n`; n=3 gives 6 coefficients including trefoil) | One beam tilt per run: averaged phase-difference image, exhaustive search over tilt and coma image shift, significance test | Per exposure group: tilt and trefoil, on by default | Beam tilt and trefoil in the forward model, from input, not fitted |
| Even aberrations | Per optics group: per-pixel 2x2 system for the phase offset, even Zernikes to n=4 (9 coefficients: Cs error, tetrafoil, group astigmatism) | None | Cs and tetrafoil, off by default | Cs, astigmatism and phase shift in the model; tetrafoil is a stub |
| Anisotropic magnification | Per optics group 2x2 matrix from a linear fit using reference gradients | Movie-level correction with user-supplied values | Optional in Global CTF Refinement | Affine resampling of the simulated image, from input |
| Ewald sphere | `relion_reconstruct --ewald` only, marked developmental: single-sideband P/Q over sectors, particle-diameter mask and weight; not in refinement | `reconstruct3d` only: double insertion with the complex CTF and its conjugate | Double insertion, about twice the reconstruction cost; both curvature signs must be tried | Intrinsic: Fresnel propagation through z-slices |
| Reference | The particle's own half map, masked, FSC-weighted, low-resolution cutoff (30 A for defocus, 20 A for aberrations) | Single reference, SSNR-weighted, 30 A low cutoff during the scan | Half maps | The volume is the unknown |

How RELION applies the result is worth copying. Odd aberrations are removed
by demodulating the particle Fourier transform before alignment and
reconstruction (`obsModel.demodulatePhase`), so the CTF stays real. Even
aberrations enter as an additive per-pixel phase inside the CTF
(`getGammaOffset`). The magnification matrix enters the projection geometry
(`applyAnisoMag`).

cryoSPARC staff report on their forum that using the Ewald-aware
image-formation model in alignment was not significant on most real
datasets; the gain is in the reconstruction.

Two points to keep in mind when cisTEM is used as a comparator.
`beam_tilt_group` exists in its star format, but `refine_ctf` does not split
by it. In `src/gui/Generate3DPanel.cpp` (lines 617-624) "apply Ewald
correction = yes" maps to `correct_ewald_sphere = 0`, which reads inverted.

Code read:

- RELION: `src/jaz/single_particle/ctf/` (`ctf_refiner.cpp`,
  `modular_ctf_optimisation.cpp`, `defocus_estimator.cpp`,
  `tilt_estimator.cpp`, `aberration_estimator.cpp`,
  `magnification_estimator.cpp`), `src/reconstructor.cpp`, `src/ctf.cpp`.
- cisTEM: `src/programs/refine_ctf/refine_ctf.cpp`, `src/core/ctf.cpp`,
  `src/core/reconstruct_3d.cpp`.
- specter: `src/specter/ghostbuster.py`, `scattering.py`, `microscope.py`,
  `imagegenerator/`.

## 3. specter

specter is a simulator plus a reconstruction (Ghostbuster) that runs gradient
descent through the simulator. The forward model rotates the real-space
volume for each particle, slices it along z, propagates the wave (projection,
first Born, Rytov, kinematic or multislice), multiplies the complex exit wave
by `exp(-i chi)` (defocus, astigmatism, Cs, beam tilt, trefoil, phase shift)
and takes `abs(psi)^2`. The reconstruction runs AdamW on the voxel values
from a zero volume with minibatches of three particles, Rytov by default,
with poses and CTF parameters read from a cryoSPARC `.cs` file.

What is new relative to the other packages: curvature, the defocus gradient
through the particle and odd aberrations are handled in one model, without
sideband splitting or a diameter-mask heuristic. The multislice mode adds
multiple scattering and the quadratic intensity term, which no Fourier-slice
method models.

What it costs: each particle needs a full volume rotation and N slice FFTs,
O(N^3) per particle against O(N^2) for a Fourier slice, with autograd memory
on top, and there is no FSC-based regularisation. The Ghostbuster paper
(arXiv 2312.08965) shows that the first-Born model is mathematically the
Ewald-sphere model. For single scattering, a curved-surface operator in
Fourier space therefore reaches the same physics at Fourier-slice cost.

## 4. When each effect matters

Estimates at 300 kV and Cs 2.7 mm. At 200 kV the tolerances are roughly 30 to
60 % tighter. 100 A is 0.01 um in SIMPLE's defocus units.

| Effect | Size that gives a pi/4 phase error |
| --- | --- |
| Defocus error | 114 A at 3 A, 79 A at 2.5 A, 51 A at 2 A, 29 A at 1.5 A |
| Beam tilt | 0.32 mrad at 3 A, 0.19 at 2.5 A, 0.10 at 2 A, 0.04 at 1.5 A |
| Anisotropic magnification | 1 % gives 0.8 rad at 3 A for an atom 75 A from the particle centre; scales with radius over resolution |

Particles sit at different heights in ice that is 300 to 500 A thick, so
per-particle defocus is the first limit.

For the Ewald sphere a simple model is used, not a measurement: over
uniformly distributed views an atom at radius r sees its height uniformly in
`[-r, r]`, and its signal is attenuated by `sinc(pi lambda r / d^2)`. The
peripheral atoms of a particle of diameter D then lose 10 % at
`d = sqrt(2 lambda D)`.

| Particle diameter | Resolution at which the periphery loses 10 % | Peripheral attenuation at 2 A |
| --- | --- | --- |
| 100 A | 2.0 A | 0.90 |
| 150 A | 2.4 A | 0.79 |
| 300 A | 3.4 A | 0.32 |
| 700 A | 5.3 A | below zero (contrast reversal) |

## 5. Where SIMPLE stands

- CTF model (`src/main/ctf/simple_ctf.f90`): `dfx`, `dfy`, `angast` and an
  additive phase. No beam tilt, no Zernike terms, no magnification matrix.
- Per-particle defocus today comes from the patch fit in `ctf_estimate`: a
  10-term polynomial surface per micrograph, evaluated at the particle
  coordinate (`doc/algorithms/ctf_estimation.md`). It follows specimen tilt
  and bending, not the height of each particle in the ice.
- Optics groups exist (`ogid`, `os_optics`, `assign_optics_groups`), but
  `os_optics` is not the analysis-time source of truth
  (`doc/policies/microscope_parameters_policy.md`).
- The Cartesian pose framework (`simple_cartesian_pose_refiner`, to become
  `cartesianft_calc` and `cartesianft_pose_opt`) provides the whitened
  residual of a particle against a fixed reference, with state-by-half
  references, and has been validated against independent oracles.
- L-BFGS-B (`src/main/opt/simple_opt_lbfgsb.f90`) is already used with
  analytic gradients in CTF estimation, the polar shift search and the
  Fourier-expanded shift search.
- `refine_motion_model` extracts per-frame particle stacks from the motion
  model; this is the start of polishing.

## 6. Plan

```text
step 1   per-particle defocus and astigmatism       <- first; waits for the Cartesian framework
step 2   beam tilt and trefoil per optics group
step 3   even aberrations, then anisotropic magnification
step 4   Ewald sphere in reconstruction only
step 5   polishing (separate track); steps 1-2 run before and again after it
```

Steps 1 to 3 share one residual: the whitened Cartesian residual at fixed
pose, against the particle's own half reference.

### 6.1 Step 1: defocus and astigmatism

With `G_w` the whitened reference projection at the particle's pose, `Y_w`
the whitened particle and `C` the CTF, the objective and its dependence on
the CTF parameters are

```text
Phi = 1/2 sum abs(C G_w - Y_w)^2
a   = abs(G_w)^2,   b = Re(conj(G_w) Y_w)        (stored once per particle)
dPhi/dp = sum (C a - b) dC/dp
```

so once `a` and `b` are stored, every objective and gradient evaluation is a
real sum over pixels with no interpolation. In SIMPLE's convention
`C = sin(chi + phi_amp)`, `chi = pi lambda s^2 (df(theta) - 1/2 lambda^2 s^2 Cs)`
and `df(theta) = 1/2 [dfx + dfy + (dfx - dfy) cos 2(theta - angast)]`, which
gives

```text
dC/d dfx    = cos(chi + phi_amp) pi lambda s^2 cos^2(theta - angast)
dC/d dfy    = cos(chi + phi_amp) pi lambda s^2 sin^2(theta - angast)
dC/d angast = cos(chi + phi_amp) pi lambda s^2 (dfx - dfy) sin 2(theta - angast)
```

These expressions were checked against central finite differences. The
optimizer is L-BFGS-B with the analytic gradient, with bounds around the
current values. The patch-surface value is the starting point.

Points to settle when the design starts, none of them decided here:

- Which parameters are per particle and which are shared per micrograph.
  RELION's default is defocus per particle and astigmatism per micrograph.
- Whether a 1D defocus scan precedes the L-BFGS-B solve. The objective is
  oscillatory in defocus for small particles; RELION offers a pre-scan and
  cisTEM uses only a scan.
- Whether to keep `(dfx, dfy, angast)` or use the linear form
  `chi ~ (dz + a1) x^2 + 2 a2 x y + (dz - a1) y^2`, which has no angle
  degeneracy at zero astigmatism.
- A prior centred on the patch-surface value, and an automatic choice
  between per-particle and per-micrograph fitting from the curvature of the
  objective. Neither exists in the other packages.
- The signed CTF is required; the phase-flipped mode needs its own handling.
- Sampling and box: unbinned particles, and a box that holds the signal
  delocalised by about `lambda defocus / d` beyond the particle on each side.

### 6.2 Step 2: beam tilt and trefoil per optics group

Two per-pixel accumulators per group, filled in the same pass as step 1, a
linear odd-Zernike fit and a significance test on the coefficients. The
correction is applied by demodulating the particle Fourier transform at
image preparation, which keeps the CTF real and serves the polar and
Cartesian branches alike. This step needs `os_optics` to become the
analysis-time source for group aberrations.

### 6.3 Step 3: even aberrations, then anisotropic magnification

Even Zernikes to n=4 enter as a per-pixel additive phase;
`ctf%eval_canonical` already takes an additive phase argument. Magnification
anisotropy changes the projection geometry in the polar sampling, the
Cartesian gather and the reconstructor, so it is deferred until a dataset
shows it.

### 6.4 Step 4: Ewald sphere in reconstruction only

Each image Fourier coefficient is a known linear combination of two volume
Fourier values, one on each cap of the sphere. That is linear in the volume,
so it fits the PCG normal equations directly and gives the least-squares
solution without a mask-and-weight heuristic. The gridding backend would use
double insertion. Alignment stays on the flat model. A by-product is the
absolute hand, from the curvature sign.

### 6.5 specter's role

Now: the independent simulator for testing steps 1 to 4. It produces
multislice stacks with known defocus spread, beam tilt and trefoil and
writes RELION star files. Later: a research track with the Singapore group
on what multislice adds beyond single scattering for large or dense
particles, started from a converged SIMPLE map and poses, not from zero.

## 7. Open decisions

- Promote `os_optics` to the analysis-time source of truth for group
  aberrations (needed for step 2).
- Ewald sphere: PCG operator first, or double insertion in the gridding
  backend first.
- Anisotropic magnification: in scope, or deferred until needed.

## 8. References

- Zivanov et al. (2018) eLife 7:e42166; Zivanov, Nakane and Scheres (2020)
  IUCrJ 7:253 (RELION CTF refinement and aberrations).
- Grant, Rohou and Grigorieff (2018) eLife 7:e35383 (cisTEM).
- Wolf, DeRosier and Grigorieff (2006) Ultramicroscopy 106:376; Russo and
  Henderson (2018) Ultramicroscopy 187:26 (Ewald sphere).
- Yeo et al., Ghostbuster, arXiv 2312.08965.
- cryoSPARC guide: Global CTF Refinement, Homogeneous Refinement; forum
  thread "Ewald sphere correction algorithm details".
