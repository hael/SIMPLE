# Ghostbuster-compatible scattering simulator

## Contract

Add a public `simple_exec prg=simulate_scattering` program that generates cryo-EM detector intensities from either a PDB file or a sampled 3D electrostatic potential. PDB input is converted with SIMPLE's existing `atoms%pdb2mrc` path. The program writes an MRC image stack and the exact orientations/CTF metadata used to generate it.

This change implements the forward models used by Ghostbuster, not its RAF phase-retrieval or DOFS reconstruction optimizers. Supported forward models are:

- `first_born`: the linear, kinematic slice sum propagated to the exit plane;
- `multislice`: nonlinear slice transmission with Fresnel propagation between slices.

Both models apply an astigmatic objective-lens phase, an optional Gaussian B-factor envelope, intensity detection, exact Poisson shot noise, and optional additive Gaussian detector noise. Defocus is measured from the object exit plane.

## Implementation plan

1. Put the wave-propagation numerics in an image-domain `scattering_simulator` type with owned complex FFT work buffers and explicit lifecycle methods.
2. Keep PDB conversion, orientation generation, volume rotation, metadata, and stack I/O in a simulation commander.
3. Expose one UI entry with exactly one of `pdbfile` and `vol1` required.
4. Register deterministic analytic tests in the existing `unit_image` test process; do not add an executable or CTest entry.
5. Document the physical conventions and command examples in the user guide.

## Review findings

- A real-to-complex image FFT cannot retain a general complex electron wave. The simulator therefore owns two aligned complex-to-complex FFTW buffers.
- The absorptive fraction belongs to the complex object potential and is not the same quantity as CTF amplitude contrast.
- SIMPLE's conventional projection simulator is not reusable for this calculation because it integrates the volume before propagation.
- The public program remains orchestration-only; the first-Born and multislice algorithms live in the image domain.

## Validation criteria

- A zero potential produces unit intensity for both models.
- A spatially constant, real potential produces the analytic first-Born intensity and unit multislice intensity when lens phase and noise are disabled.
- A spatially constant absorptive potential produces the analytic exponential attenuation under multislice.
- The parameter/UI registries and Fortran source descriptors pass their static checks.
- The focused `unit_image` sub-suite, incremental build, install, and fast gate pass.
- A small PDB-to-stack smoke test completes and writes both the image stack and orientation table.

## Status

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
