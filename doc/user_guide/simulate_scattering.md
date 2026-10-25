# simulate_scattering — Simulate Wave-Scattering Images

`simple_exec prg=simulate_scattering` generates boxed cryo-EM detector intensities from a sampled
three-dimensional electrostatic potential, or from atomic coordinates converted to one, with the wave models
of `reconstruct3D_wave`. It is a data simulator: the images follow Kirkland's sign conventions (a positive
potential is dark at underfocus) and the defocus refers to the specimen midplane, the centre of the box, as in
SIMPLE's CTF fits and in specter.

## Models

`wave_model` selects the image formation. With `sigma` the interaction constant, `a` the pixel size and
`alpha = fraca` the amplitude contrast, the potential enters as `c V` with `c = i sigma a (sqrt(1 - alpha^2) + i alpha)`.

- `projection`: the exit wave `exp(c u)` of the projected potential `u`, sampled on the central section.
- `first_born`: `1 + c u`, with `u` the volume's transform on the Ewald paraboloid of the sign `ewald`
  selects (`neg` samples `kz = -lambda |k|^2 / 2`, specter's negative sign; `pos` the opposite one; `no` the
  central section).
- `rytov`: `exp(c u)` on the paraboloid; the exponential is taken where the wave leaves the specimen, half a
  box downstream of the midplane, and the wave is carried back to the midplane for the lens.
- `multislice` (default): slice-by-slice transmission `exp(c v_z)` and Fresnel propagation
  `exp(-i pi lambda d |k|^2)` through the rotated potential; `ewald=neg` sends the beam in through the first
  section, `pos` through the last. It retains multiple scattering.

All four models run on the operator of `reconstruct3D_wave`: the PCG lattice's Kaiser-Bessel sampling of the
volume's transform, deapodised, and the image formed on a grid twice as fine. The single-scatter models sample a
section of the transform on the native Nyquist disk; multislice gathers the particle's rotated cube from the
same lattice (its native spectrum on the Nyquist sphere), slices it along the beam, and runs the recursion on the
fine grid, so its rotation is the exact transpose of the reconstruction's and needs no real-space resampling.
Particle shifts and phase shifts apply to every model. The conventions, with their comparison against specter,
are in `doc/policies/3D/wave_model_policy.md`; [reconstruct3D_wave](reconstruct3D_wave.md) reconstructs from such
images with the same models.

After the specimen the wave passes the objective lens of SIMPLE's CTF: the phase of `defocus`, the
astigmatism, `cs` and the phase shift (on the scattered waves only), with a Gaussian amplitude envelope
`exp(-bfac |k|^2 / 4)`. The recorded image is the relative intensity `|psi|^2`, vacuum 1.

When `total_dose` is positive, the intensity is converted to expected electrons per pixel using
`total_dose*smpd^2`, sampled from an exact Poisson distribution, and converted back to intensity units.
`noise_sdev` adds independent Gaussian detector noise after dose normalization.

## PDB input

SIMPLE's `pdb2mrc` implementation converts and centers the atomic model before scattering simulation:

```sh
simple_exec prg=simulate_scattering \
  pdbfile=molecule.pdb smpd=1.3 vol_dim=192 \
  nptcls=100 wave_model=multislice ewald=neg \
  kv=200 cs=2.0 fraca=0.1 defocus=1.5 dferr=1.5 astigerr=0.2 \
  bfac=64 total_dose=20 \
  outstk=scattering_images.mrcs outfile=scattering_oris.txt
```

The generated potential is `scattering_potential.mrc` unless `outvol` is supplied. `vol_dim` is a cubic pixel
dimension. It must be large enough to contain the molecule.

## Potential-volume input

To reuse an existing electrostatic potential:

```sh
simple_exec prg=simulate_scattering \
  vol1=potential.mrc smpd=1.3 nptcls=100 \
  wave_model=rytov ewald=neg total_dose=0 noise_sdev=0 \
  outstk=rytov_images.mrcs outfile=rytov_oris.txt
```

The input must be a cubic MRC volume whose values are electrostatic potential samples in volts. Supply exactly
one of `pdbfile` and `vol1`.

## Orientations and microscope settings

Without `oritab`, the program draws random orientations and CTF values from `defocus`, `dferr`, and
`astigerr`. With `oritab`, it uses the supplied orientations, shifts and any `dfx`, `dfy`, `angast` and
`phshift` fields. If the defocus fields are absent, they are generated. The output orientation table always
records the truth used for every image. `fraca`, `kv` and `cs` are the program's keys.

The defaults follow the Ghostbuster simulation regime where practical: `bfac=64` and `total_dose=20`. Voltage
and amplitude contrast retain SIMPLE's normal defaults unless `kv` and `fraca` are given.

## Practical limits

Multislice gathers the particle's rotated cube from the lattice and then performs two-dimensional complex FFTs on
the twice-finer grid for every slice and propagation step; it is considerably more expensive than the
single-scatter models, which cost a few such transforms per particle. Start with a small box and a few
particles when checking a setup.
