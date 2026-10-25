# simulate_particles — Simulate Boxed Particle Images

`simple_exec prg=simulate_particles` projects a molecular density volume to
generate boxed particle images with simulated shifts, optional CTF effects,
and noise. It writes both the image stack and its orientation/CTF metadata.

## Inputs and setup

Use an installation configured as described in
[Installation](../installation.md), with `simple_exec` on `PATH`.
Supply a centered cubic MRC or SPIDER density volume; this command does not
read PDB coordinates directly. Set `smpd` to the volume's actual sampling in
Å/pixel. Changing this option does not resample the volume.

The examples assume a 192 × 192 × 192 MRC volume at **1.3 Å/pixel**, with
the molecule fitting within a 180 Å mask. Replace
`/absolute/path/reference.mrc` with your volume and adjust these settings
for your molecule. Run each example in its own fresh directory.

## Examples

### Noisy particles with CTF variation

```shell
mkdir particles_ctf
cd particles_ctf

simple_exec prg=simulate_particles \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 mskdiam=180 pgrp=c1 nptcls=200 snr=0.5 sherr=2 \
  ctf=yes kv=300 cs=2.7 fraca=0.1 \
  defocus=2 dferr=0.2 astigerr=0.05 bfac=200 bfacerr=0 \
  outstk=particles.mrcs outfile=particle_truth.txt mkdir=no nthr=1
```

This generates 200 images at random orientations, with shifts bounded by
2 pixels per image axis. The defocus components and astigmatism vary per
particle; the chosen values are saved with each orientation record.

### A mixture of states

```shell
mkdir particles_states
cd particles_states

simple_exec prg=simulate_particles \
  vol1=/absolute/path/state1.mrc vol2=/absolute/path/state2.mrc vol3=/absolute/path/state3.mrc \
  nstates=3 state_pops=0.5,0.3,0.2 \
  smpd=1.3 mskdiam=180 pgrp=c1 nptcls=1000 snr=0.5 sherr=2 \
  ctf=yes kv=300 cs=2.7 fraca=0.1 defocus=2 dferr=0.2 astigerr=0.05 \
  outstk=particles.mrcs outfile=particle_truth.txt mkdir=no nthr=1
```

This generates 1000 images of three conformations or compositions of one
molecule: 500 of `state1.mrc`, 300 of `state2.mrc` and 200 of `state3.mrc`.
Each image is projected from the volume of its particle's state, and the
state is written with the particle's orientation, shift and CTF parameters,
so the output orientation document is the truth of the mixture. The volumes
must share box, sampling and frame.

### Noisy particles without CTF

```shell
mkdir particles_no_ctf
cd particles_no_ctf

simple_exec prg=simulate_particles \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 mskdiam=180 pgrp=c1 nptcls=200 snr=0.5 sherr=2 \
  ctf=no bfac=0 bfacerr=0 \
  outstk=particles_no_ctf.mrcs outfile=particle_truth.txt mkdir=no nthr=1
```

`ctf=no` bypasses CTF filtering in this command. A B-factor is still applied
when supplied, independently of the CTF switch, so this example explicitly
sets `bfac=0` as well. Gaussian noise remains enabled by `snr=0.5`.

### Particles under a wave model

```shell
mkdir particles_wave
cd particles_wave

simple_exec prg=simulate_particles \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 mskdiam=180 pgrp=c1 nptcls=200 snr=10 sherr=2 \
  ctf=yes kv=300 cs=2.7 fraca=0.1 defocus=1 dferr=0.2 astigerr=0.05 \
  wave_model=rytov total_dose=40 \
  outstk=particles.mrcs outfile=particle_truth.txt mkdir=no nthr=8
```

With `wave_model` the images are the wave model's instead of the Fourier
projection's (see [Wave mode](#wave-mode)). `snr=10` adds no detector noise,
so the only noise is the shot noise of 40 electrons per square Angstrom.

## Controls and defaults

These defaults describe the current standalone simulation path; several
commander defaults differ from the generic parameter or help-text values.

| Option | Meaning | Effective default or requirement |
| --- | --- | --- |
| `vol1` | Input density volume; with `nstates=K`, `vol1` to `volK`, one per state | Required |
| `nstates` | Number of states, each with its own volume | 1 |
| `state_pops` | Comma-separated population fractions of the states, summing to 1 | Equal populations |
| `smpd` | Sampling in Å/pixel | Supply the actual input sampling; standalone parameter initialization is 1.3 |
| `nptcls` | Number of output particle images | Required |
| `snr` | Signal-to-noise setting | Required; use a positive value |
| `mskdiam` | Soft spherical mask diameter in Å | Box-derived mask when omitted |
| `even` | Generate quasi-even projection directions | `no`: random orientations unless `oritab` is supplied |
| `oritab` | SIMPLE orientation table | Omitted: generate orientations |
| `sherr` | Uniform shift half-width in pixels for random orientations | 0 |
| `trs` | Uniform shift half-width in pixels with `even=yes` | 0 |
| `pgrp` | Point-group setting | `c1` |
| `ctf` | Apply CTF filtering | `yes`; use `no` to bypass it (no other value is accepted) |
| `kv` | Acceleration voltage in kV | 300 |
| `cs` | Spherical aberration in mm | 2.7 |
| `fraca` | Amplitude contrast fraction | 0.1 |
| `defocus` | Nominal underfocus in µm | 2 |
| `dferr` | Variation of `dfx` around nominal defocus, in µm | 1.5 |
| `astigerr` | Variation of `dfy` around `dfx`, in µm | 0.5 |
| `deftab` | Existing SIMPLE CTF-parameter table | Omitted: generate CTF parameters when needed |
| `bfac` | Image B-factor in Å² | No B-factor filtering unless this key is supplied |
| `bfacerr` | Uniform B-factor variation half-width in Å² | 0; used only when `bfac` is supplied |
| `outstk` | Particle-stack filename | `simulated_particles.mrc` for MRC input |
| `outfile` | Orientation/CTF metadata filename | `simulated_oris.txt` |
| `mkdir` | Create a numbered execution directory | `no` |
| `nthr` | OpenMP thread count | Runtime OpenMP maximum when omitted; examples use 1 |
| `wave_model` | Wave mode: `projection`, `first_born`, `rytov` or `multislice` | Omitted: the Fourier-projection simulator |
| `ewald` | Ewald-sphere curvature of the wave model: `neg`, `pos` or `no` | `neg`, `no` for `projection`; refused without `wave_model` |
| `total_dose` | Electrons per square Angstrom: the scale link and the shot noise | Required with `wave_model`; refused without it |

### Orientations, shifts, and CTF parameters

With neither `even=yes` nor `oritab`, the command generates random Euler
orientations and shifts using `sherr`. With `even=yes`, it instead generates
quasi-even projection directions, random in-plane angles, and shifts using
`trs`. Set `trs` explicitly if you want shifted quasi-even projections.

To reuse existing poses, supply a SIMPLE orientation table with `oritab=`,
set `nptcls` to its record count, and retain `even=no`. The table's shifts
are used as supplied; `sherr` does not perturb them. `even=yes` regenerates
poses even if an orientation table is also supplied.

When CTF metadata is already present, it is retained instead of randomly
generated. `deftab` can supply CTF records; `oritab` is read afterward and
overrides those records when both are supplied. Otherwise, with `ctf=yes`,
the generator draws positive `dfx` values within `defocus ± dferr`, then
positive `dfy` values within `dfx ± astigerr`, with random astigmatism angles.

### States

With `nstates=K` greater than 1, every particle has a state from 1 to K and is
projected from that state's volume. If `oritab` gives the particles a state
(a `state=` field; a project file always does), those states are used, a
state outside 1 to K is an error, and so is `state_pops` beside them (the
populations would not be applied). Otherwise (no `oritab`, or one without
states) the states are drawn with the fractions of `state_pops`: one fraction
per state, none negative, summing to 1 within 0.001, or else an error. They are
drawn in exact counts: each fraction times `nptcls`, rounded, with what rounding leaves
over given to (or taken from) the states in order of decreasing fraction. The
order of the states over the particles is then shuffled with the seeded random
generator, after the orientations and CTF parameters are drawn. The K padded
volumes are held in memory together, so the memory use grows with K. A state
that ends with no particle (a population too small for `nptcls`) is an error. With
`nstates` absent or 1, every particle has state 1 and the output is the same as
with one volume.

The command always soft-masks each input volume before projection. If
`mskdiam` is omitted, the default diameter is `(box - 6) * smpd` Å before
conversion to the supported pixel radius. Invalid or oversized mask diameters
fall back to the box-derived mask with a warning.

### Wave mode

`wave_model` switches the image formation to SIMPLE's wave operator, the one
`reconstruct3D_wave` and `simulate_scattering` use. Each state's volume, masked
as above, is read as an ordinary SIMPLE map and converted to a potential in
volts by the scale link `kappa = 2 g sigma a N` (`g = sqrt(total_dose) a` the
gain, `sigma` the interaction constant at `kv`, `a` the pixel size, `N` the
box). For every particle the operator predicts the expected electron counts
`E = total_dose a^2 I` (`I` the relative intensity) at its pose, shift, CTF and
lens B-factor (`bfac`, `bfacerr` act on the lens envelope of each particle).
Each count is drawn from a Poisson distribution of mean `E` (shot noise); for
`snr` below 5, Gaussian detector noise of standard deviation
`sqrt(var(E) / snr)` follows, the same rule as the Fourier-projection mode.
The particle is written as `x = -(counts - D a^2) / sqrt(D a^2)` (`D` the dose): contrast
inverted and normalised to unit shot noise, the units `reconstruct3D` expects.
A stack made this way and reconstructed with `reconstruct3D` (PCG) returns its
map's scale to within a few percent at low and middle resolution, as long as
the particles' background is free of signal: a defocus whose delocalisation
reaches past the mask into the background normalisation lowers the scale.

Wave mode needs `ctf=yes` (a phase object in focus has almost no contrast) and
`total_dose`; `projection` refuses a curvature, `multislice` needs one. Because
the scale link fixes the map's scale in volts, the dose does not change the
particles' signal-to-noise ratio, only how strongly nonlinear the image is.

The two modes differ in scale: the Fourier-projection mode projects the volume
padded to twice the box and normalises over the padded box, without
deapodisation, so its noiseless particle is about half the wave mode's weak
limit. That mode's absolute scale is not meant to be used: its noise is set
relative to the image, and reconstruction normalises particles to unit noise.

Cost on a 24-core machine at box 128, 2000 particles: about 2.5 s for
`projection`, `first_born` and `rytov` (1.2 to 1.4 ms per particle), about
95 s for `multislice` (47 ms per particle), against about 7 s for the
Fourier-projection mode.

## Outputs

| File | Contents |
| --- | --- |
| `outstk` (default `simulated_particles.mrc` for MRC input) | One boxed projection image per particle |
| `outfile` (default `simulated_oris.txt`) | SIMPLE records containing generated or supplied orientations, shifts, states, and CTF parameters when applicable |

The output image box follows the input volume box. When `outstk` is omitted,
the default filename extension follows the input volume format. An explicit
`outstk` selects the output format through its extension: `.mrc` or `.mrcs`
for MRC, and `.spi` for SPIDER. Use separate working directories or unique
output names to preserve results from repeated runs.

## Limitations

Without `wave_model` this is a simplified Fourier-projection simulator, not a
multislice electron-scattering model for heavy-atom specimens; with it, see
[Wave mode](#wave-mode). It generates boxed
images rather than full micrographs or local specimen motion.
The current path writes one image per particle; the advertised `nframes`
option does not generate a sequence of particle movie frames here.

The implementation adds Gaussian detector noise when `snr < 5`; values
of 5 or higher skip noise addition. CTF and B-factor filtering precede this
noise. Do not interpret the older help-text description of a 20%/80% noise
split as the behavior of this particle generator.

For clean projections to place in a micrograph, use [reproject](reproject.md).
For a full-frame movie and aligned average, use
[simulate_movie](simulate_movie.md). Built-in fixtures and quantitative checks
are documented in [Simulated validation workflows](simulated_workflows.md).
