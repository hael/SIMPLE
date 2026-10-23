# simulate_scattering — Simulate Wave-Scattering Images

`simple_exec prg=simulate_scattering` generates boxed cryo-EM detector intensities with the forward physics used to create Ghostbuster test data. It accepts either atomic coordinates or a sampled three-dimensional electrostatic potential. This is a data simulator; it does not run Ghostbuster's RAF phase retrieval or DOFS reconstruction.

## Models

`scattering_model=first_born` implements the kinematic slice sum

\[
\psi_\mathrm{exit}\simeq 1+i\sigma_e\sum_m v_m\otimes p_m,
\]

where every potential slice is propagated from its depth to the object exit plane. `scattering_model=multislice` alternates complex slice transmission and Fresnel propagation, retaining multiple-scattering effects.

The object potential is

\[
v=\left(\sqrt{1-\alpha^2}+i\alpha\right)V,
\]

where `absorption` supplies \(\alpha\). After the object, the simulator applies an astigmatic objective-lens phase with `defocus`, `cs`, and a Gaussian `bfac` amplitude envelope. Defocus is measured from the object exit plane. The recorded image is \(|\psi|^2\).

When `total_dose` is positive, the intensity is converted to expected electrons per pixel using `total_dose*smpd^2`, sampled from an exact Poisson distribution, and converted back to intensity units. `noise_sdev` adds independent Gaussian detector noise after dose normalization.

## PDB input

SIMPLE's `pdb2mrc` implementation converts and centers the atomic model before scattering simulation:

```sh
simple_exec prg=simulate_scattering \
  pdbfile=molecule.pdb smpd=1.3 vol_dim=192 \
  nptcls=100 scattering_model=multislice \
  kv=200 cs=2.0 defocus=1.5 dferr=1.5 astigerr=0.2 \
  absorption=0.1 bfac=64 total_dose=20 \
  outstk=scattering_images.mrcs outfile=scattering_oris.txt
```

The generated potential is `scattering_potential.mrc` unless `outvol` is supplied. `vol_dim` is a cubic pixel dimension. It must be large enough to contain the molecule.

## Potential-volume input

To reuse an existing electrostatic potential:

```sh
simple_exec prg=simulate_scattering \
  vol1=potential.mrc smpd=1.3 nptcls=100 \
  scattering_model=first_born total_dose=0 noise_sdev=0 \
  outstk=born_images.mrcs outfile=born_oris.txt
```

The input must be a cubic MRC volume whose values are electrostatic potential samples in volts. Supply exactly one of `pdbfile` and `vol1`.

## Orientations and microscope settings

Without `oritab`, the program draws random orientations and CTF values from `defocus`, `dferr`, and `astigerr`. With `oritab`, it uses the supplied orientations and any `dfx`, `dfy`, and `angast` fields. If those CTF fields are absent, they are generated. The output orientation table always records the truth used for every image.

The defaults follow the Ghostbuster simulation regime where practical: `absorption=0.1`, `bfac=64`, and `total_dose=20`. Voltage retains SIMPLE's normal default unless `kv` is given.

## Practical limits

The implementation rotates the full potential once per particle and then performs a two-dimensional complex FFT for every slice and propagation step. Multislice is therefore considerably more expensive than projection simulation. Start with a small box and a few particles when checking a setup.
