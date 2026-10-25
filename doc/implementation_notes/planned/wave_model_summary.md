# Wave-model image formation in SIMPLE: where we are and what comes next

Status: summary for review, 2026-10-24. Nothing below is started unless it says
so.

## Terms

- **Wave model**: predicting a cryo-EM image from the specimen's electrostatic
  potential with wave optics, as specter's Ghostbuster does. The beam picks up
  a phase in the specimen, passes the objective lens, and the detector records
  the squared modulus of the wave. Four variants: projection (thin specimen),
  first Born and Rytov (single scattering, with the Ewald-sphere curvature),
  and multislice (multiple scattering, slice by slice).
- **Flat model**: SIMPLE's existing image model, CTF times the central section
  of the map. It is the weak-signal limit of every wave model.
- **Ewald-sphere curvature**: an image samples the map's transform on a curved
  surface, not a plane. It limits large particles, low voltages and high
  resolution. Its sign is a property of the microscope and has to be tried.
- **Scale link**: every map SIMPLE writes is the specimen's potential in volts
  times one number, `kappa = 2 g sigma a N` (`g` the square root of the dose
  per pixel, `sigma` the interaction constant at the voltage, `a` the pixel
  size, `N` the box). It follows from the flat model being the weak-signal limit
  of the wave model, and lets any SIMPLE map serve as a potential.
- **specter**: the Python package of our collaborators that contains
  Ghostbuster; used here as an independent reference implementation.

## 1. Done: the wave operators and reconstruct3D_wave

An unattended run on the Dell (2026-10-23/24) added the wave model to SIMPLE
as Fortran, built on SIMPLE's own machinery rather than translated from
specter: the PCG reconstructor's interpolation lattice (now with the Ewald
curvature), the CTF code, the particle reading of reconstruction, the PCG
start map, density and prior. 28 files changed, 17 of them new. The changes
are in `~/src/SIMPLE-for-codex`, not yet in the repository.

What it adds:

- the forward model, its derivative and its transpose (the "backprojection")
  for all four variants, with either curvature sign;
- `reconstruct3D_wave`, a standalone shared-memory program that reconstructs
  half maps from a project with known poses and CTF: per half, the flat PCG map
  as the start, then a nonlinear conjugate-gradient solve under the wave model
  with SIMPLE's FSC-based prior;
- `simulate_scattering` on the same operator, with its sign and defocus errors
  corrected (it imaged the negative of the potential and measured defocus from
  the exit face).

How it was checked: adjoint tests (agreement 2e-8 to 2e-6), gradients against
finite differences (1e-4 or better), the weak-signal limit against SIMPLE's
existing projector (2e-6), identical results at 1, 11, 16 and 24 threads, PCG
results byte-identical to before, and agreement with specter's images to about
1% (the interpolation error of SIMPLE's three-point window).

Ground truth, datasets simulated by specter with multislice:

| Dataset | Band | Flat PCG | Wave model, correct sign | Wrong sign | specter |
| --- | --- | --- | --- | --- | --- |
| RuBisCO, 300 kV, 6 000 particles | 3-2 A | 0.730 | 0.762 (Rytov) | 0.707 | 0.765 |
| Ribosome, 100 kV, 4 000 particles | 3-2.4 A | 0.650 | 0.810 (Rytov), 0.779 (first Born) | 0.21 | 0.833 |
| Pt nanoparticle, 2 000 particles | atoms | all 469 found, 0.035 A error | all found, 0.029 A error | | |

(Fourier shell correlation against the true map in the band shown.) On the
ribosome, `reconstruct3D_wave` took 9.5 minutes against 1.5 for PCG on 24
threads; one pass over the particles costs 1.7 to 3 times a PCG accumulation.

Settled after the run:

- SIMPLE's noise normalisation stays. It overestimates the noise when the
  CTF spreads signal into the background ring (up to 10% at 100 kV), which
  makes the absolute scale approximate but leaves the map unchanged.
- The ribosome's last band stays 0.02 below specter's; it appears only where
  the CTF spreads signal beyond the box margin.
- The fast-gate failure seen during the run (`unit_pftc_align2D3D`) is fixed on
  master.
- `simulate_scattering` stays.

Open from the run, not yet decided:

- The prior, built from the flat start maps, holds back the band near the
  curvature limit; a ten times weaker prior recovered most of it.
- `reconstruct3D_wave` inherits `objfun=euclid` and stops on projects without
  noise spectra unless `objfun=cc` is given.
- A seven-point interpolation window cut the 1% error to 5e-5; an option for
  later.
- `rotvol`, SIMPLE's volume rotation, does not correct the interpolation
  roll-off and attenuates density away from the centre by 8 to 10%. The wave
  code no longer uses it; other callers may want a look.

## 2. Ready to launch: STAR origins

Found by the same run: every STAR import and export path copies SIMPLE's
particle shifts (pixels) into RELION's Angstrom origins without the pixel size
(1.5 pixels off at 1.3 A per pixel), `export_relion` writes no poses, and a
STAR file without an optics block fails with an unhelpful message. A two-phase
Dell run fixes all paths, makes `export_relion` write poses, refuses the old
RELION 3.0 layout with a clear message, and checks positions both ways through
specter. Plan: `doc/implementation_notes/planned/star_origins_angstrom.md`.
It runs as the first phase of the combined run `star_wave` (see "Order of
work").

## 3. Proposed next run: simulate_particles, speed, clean-up

**A wave option in `simulate_particles`.** `wave_model=` (and `ewald=`) on the
existing program; without the key it behaves exactly as today. In wave mode:

- the input is a SIMPLE map, turned into a potential with the scale link from
  `total_dose`, the voltage, the pixel size and the box. Any reconstruction can
  then be imaged with the wave model; `simulate_scattering` stays the route from
  a PDB model;
- noise in two steps, as a counting camera records it: shot noise, Poisson at
  `total_dose` (the dose that also sets the scale) on the expected electron
  counts; then detector noise, Gaussian in real space, with `snr` keeping its
  meaning (the noiseless image's variance over the detector noise's, as
  `add_gauran` sets it today, none for `snr` of 5 or more). The particles
  come out in SIMPLE's convention (contrast inverted, normalised);
- multiple states, `recthres`, parallel parts and the B-factor options work as
  today, since only the per-particle image call changes;
- the scale link gets one home (a function), used by `reconstruct3D_wave` and
  `simulate_particles`.

**Speed** (measured on the Dell before and after):

- the lens transfer function is rebuilt for every particle on every pass with
  several transcendental calls per fine-grid point; the parts that do not
  change between particles are computed once, leaving one sine and cosine per
  point;
- the intensity and its derivative are real fields transformed with complex
  FFTs; real-to-complex transforms halve those;
- work arrays allocated per particle inside the parallel loops become
  per-thread buffers allocated once;
- the shift phase is computed per pixel; separable phase vectors do it
  cheaper.

**Clean-up to SIMPLE's rules**: shorter code comments (the type headers run to
20 lines against the limit of five), no umbrella-module imports in the new
modules, the image-writing loop in `simulate_scattering`'s commander moved into
the operator, and code comments pointing at the policy
(`doc/policies/3D/wave_model_policy.md`) instead of the plan.

Decided for this run (2026-10-24): Poisson shot noise at the dose followed by
Gaussian detector noise in real space; the input of wave mode is an ordinary
SIMPLE map converted with the scale link. The run's plan is
`doc/implementation_notes/planned/wave_simulate_speed.md`.

## 4. Later: coupling to refinement (to plan with Ruben)

No PDB model is needed: through the scale link the map of the previous
iteration is the potential. Where the wave model enters:

- **Reconstruction**: `rec_backend=wave` in `reconstruct3D` and `refine3D`,
  with distributed execution (workers publish gradient sums, the master steps).
- **Continuous pose refinement** (`refine=cont`) and **CTF refinement**:
  per-particle and gradient-based, where the full wave model fits. They need
  the per-particle forward model (it exists) and its derivatives with respect
  to pose and defocus (the interpolation window's gradient exists).
- **Discrete orientation search**: wave images depend on each particle's
  defocus, so they cannot be precomputed as references. To first order an image
  is two curved slices of the map combined with per-particle lens factors; the
  search would need two slices per orientation instead of one. At 300 kV on
  ordinary proteins the two coincide; the search can stay flat at first.
- Not yet scheduled: tomography, beam tilt and trefoil in the lens, an
  Ewald-aware gridding backend.

## Order of work

The work stays in `~/src/SIMPLE-for-codex` until it goes to the repository.

1. The `wave_operators` changes are there, committed on `811934ca5` (2026-10-24).
2. One unattended run, `star_wave`, carries out section 2 and then section 3
   (run plan `star_origins_and_wave_simulate_run.md`).
3. Plan the refinement coupling (section 4) with Ruben.
