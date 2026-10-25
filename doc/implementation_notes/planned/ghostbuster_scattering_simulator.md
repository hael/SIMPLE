# Ghostbuster-compatible scattering simulator

## Update, 2026-10-24: conventions corrected, simulator on the shared operator

`simulate_scattering` now runs on the wave operators that `reconstruct3D_wave` reconstructs with
(`doc/implementation_notes/planned/ghostbuster_operators.md`, its report beside it, and
`doc/policies/3D/wave_model_policy.md`). What changed for this program, and why:

- **Signs.** The first implementation conjugated both the Fresnel propagator (`exp(+i pi lambda d k^2)`) and
  the lens (`exp(+i chi)`) relative to Kirkland, while keeping the transmission `exp(+i sigma V)`. Its images
  were therefore those of the negated potential: a positive potential came out bright at underfocus, and the
  multiple-scattering terms were those of `-V`. Measured against specter on the same potential and pose, the
  simulator's image of `V` correlated with specter's image of `-V` at 0.98 to 1.000 and with specter's image of
  `V` at -0.90 to -0.93. Only the phase was inverted: the absorption had the right sign. The program now uses
  Kirkland's signs throughout, as specter does.
- **Defocus reference.** Defocus was measured from the exit plane; it now refers to the specimen midplane (the
  box centre), as SIMPLE's CTF fits and specter do. The two differ by half the box depth (48 A in the test case).
- **Keys.** `scattering_model` is now `wave_model` (`first_born` and `multislice` keep their names;
  `projection` and `rytov` are new), with `ewald=neg|pos|no` for the curvature sign (for multislice, the beam
  direction). `absorption` is retired: it is the same quantity as the amplitude contrast, so `fraca` (the same
  default, 0.1) sets it. With the object potential `(sqrt(1 - alpha^2) + i alpha) V`, the weak-phase limit of
  every model is SIMPLE's CTF with amplitude contrast `alpha` (tested to `1e-5`), which revises the review finding
  of 2026-10-23 below that the two differ. No aliases.
- **Rotation.** The volume is no longer rotated in real space (`rotvol`) per particle. Every model samples the
  potential's transform on the padded lattice of SIMPLE's PCG reconstruction with its Kaiser-Bessel window,
  deapodised; multislice gathers the particle's rotated spectrum cube from the same lattice. The data are made by
  the exact operator the reconstruction inverts. (`rotvol` did not correct the window's roll-off: a weak
  potential's image came out 8 to 10% weaker off the box centre.)
- **Lens.** SIMPLE's CTF phase (`ctf%eval_phase`, the same astigmatism and Cs convention as every CTF in SIMPLE),
  with the phase shift acting on the scattered waves only, and the B-factor envelope `exp(-bfac |k|^2 / 4)`.
- **Shifts and phase shifts** of the orientation table are now applied, for every model, and the output table
  records the CTF parameters used.
- **Code.** `simple_scattering_simulator.f90` and its tester are deleted. The FFT buffers became `wave_field`
  (`src/main/image/simple_wave_field.f90`); the lens, the intensity and the dose and detector noise
  (`detect_intensity`) live in `wave_imaging` (`src/main/image/simple_wave_imaging.f90`), multislice in
  `wave_multislice` (`src/main/volume/simple_wave_multislice.f90`), the particle loop in `wave_projector`
  (`src/main/volume/simple_wave_projector.f90`); the Poisson sampler is `poisson_rnd` in `simple_rnd`, the
  wavelength and interaction constant are in `simple_math_ctf`. The FFTW plans are made under the unnamed
  critical section that every FFTW planner call of SIMPLE takes (the simulator's named section was a second lock).
- **Tests.** The multislice limits (zero, constant real and constant absorptive potentials) and the detector
  checks are kept; the first-Born constant-potential check gave way to the weak-phase identity (the linear term
  of every model is SIMPLE's CTF times the projection) and to first Born being the linearisation of Rytov, and
  the vacuum is tested for every model. New tests can see signs and the defocus reference, which the old ones could not: a weak
  sheet of potential in the midplane images as SIMPLE's weak-phase model (CTF times the projection); the same
  sheet 16 A downstream images at `df + 16 A` and not at `df - 16 A`; with the beam reversed at `df - 16 A`
  (`wave multislice` in `unit_reconstruction`). `scattering simulator` left `unit_image`.

Unchanged: the PDB route (`pdbfile`, `vol_dim`, `outvol`, `pdb2mrc`), the potential-volume route, the Poisson
dose (`total_dose`), the Gaussian detector noise (`noise_sdev`), the defaults `bfac=64` and `total_dose=20`,
the output stack and orientation table. The user guide is `doc/user_guide/simulate_scattering.md`.

## Contract

`simple_exec prg=simulate_scattering` generates cryo-EM detector intensities from either a PDB file or a sampled
3D electrostatic potential. PDB input is converted with SIMPLE's existing `atoms%pdb2mrc` path. The program writes
an MRC image stack and the exact orientations and CTF metadata used to generate it.

It implements the forward models used by Ghostbuster, not its RAF phase-retrieval or DOFS reconstruction
optimizers: projection, first Born and Rytov on the Ewald paraboloid, and multislice. All models apply the
objective lens of SIMPLE's CTF, an optional Gaussian B-factor envelope, intensity detection, exact Poisson shot
noise, and optional additive Gaussian detector noise. Defocus refers to the specimen midplane.

## First implementation (2026-10-23), as recorded then

The sections below record the first implementation. Its conventions (exit-plane defocus, the conjugated
propagator and lens, the `absorption` key) were corrected on 2026-10-24, above.

### Implementation plan

1. Put the wave-propagation numerics in an image-domain `scattering_simulator` type with owned complex FFT work buffers and explicit lifecycle methods.
2. Keep PDB conversion, orientation generation, volume rotation, metadata, and stack I/O in a simulation commander.
3. Expose one UI entry with exactly one of `pdbfile` and `vol1` required.
4. Register deterministic analytic tests in the existing `unit_image` test process; do not add an executable or CTest entry.
5. Document the physical conventions and command examples in the user guide.

### Review findings

- A real-to-complex image FFT cannot retain a general complex electron wave. The simulator therefore owns two aligned complex-to-complex FFTW buffers.
- The absorptive fraction belongs to the complex object potential and is not the same quantity as CTF amplitude contrast.
- SIMPLE's conventional projection simulator is not reusable for this calculation because it integrates the volume before propagation.
- The public program remains orchestration-only; the first-Born and multislice algorithms live in the image domain.

### Validation criteria

- A zero potential produces unit intensity for both models.
- A spatially constant, real potential produces the analytic first-Born intensity and unit multislice intensity when lens phase and noise are disabled.
- A spatially constant absorptive potential produces the analytic exponential attenuation under multislice.
- The parameter/UI registries and Fortran source descriptors pass their static checks.
- The focused `unit_image` sub-suite, incremental build, install, and fast gate pass.
- A small PDB-to-stack smoke test completes and writes both the image stack and orientation table.

### Status

Implementation and local macOS validation are complete.

- The incremental Release build and install pass on the current `master`.
- The 15-entry fast gate passes.
- The focused scattering suite passes all 8 analytic assertions.
- An H2O PDB smoke test passes through PDB conversion, first-Born simulation,
  multislice simulation, objective-lens propagation, Poisson dose, and
  additive Gaussian noise.
- The test registry, FLEX DAG, source descriptors, and diff checks pass.

Linux/BOX validation remains outstanding. MRC stacks written incrementally can
retain a zero header minimum even when a direct scan finds all noiseless pixels
positive; the stored image data are unaffected.
