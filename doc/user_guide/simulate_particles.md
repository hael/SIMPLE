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

## Controls and defaults

These defaults describe the current standalone simulation path; several
commander defaults differ from the generic parameter or help-text values.

| Option | Meaning | Effective default or requirement |
| --- | --- | --- |
| `vol1` | Input density volume | Required |
| `smpd` | Sampling in Å/pixel | Supply the actual input sampling; standalone parameter initialization is 1.3 |
| `nptcls` | Number of output particle images | Required |
| `snr` | Signal-to-noise setting | Required; use a positive value |
| `mskdiam` | Soft spherical mask diameter in Å | Box-derived mask when omitted |
| `even` | Generate quasi-even projection directions | `no`: random orientations unless `oritab` is supplied |
| `oritab` | SIMPLE orientation table | Omitted: generate orientations |
| `sherr` | Uniform shift half-width in pixels for random orientations | 0 |
| `trs` | Uniform shift half-width in pixels with `even=yes` | 0 |
| `pgrp` | Point-group setting | `c1` |
| `ctf` | Apply CTF filtering | `yes`; use `no` to bypass it |
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

The command always soft-masks the input volume before projection. If
`mskdiam` is omitted, the default diameter is `(box - 6) * smpd` Å before
conversion to the supported pixel radius. Invalid or oversized mask diameters
fall back to the box-derived mask with a warning.

## Outputs

| File | Contents |
| --- | --- |
| `outstk` (default `simulated_particles.mrc` for MRC input) | One boxed projection image per particle |
| `outfile` (default `simulated_oris.txt`) | SIMPLE records containing generated or supplied orientations, shifts, and CTF parameters when applicable |

The output image box follows the input volume box. When `outstk` is omitted,
the default filename extension follows the input volume format. An explicit
`outstk` selects the output format through its extension: `.mrc` or `.mrcs`
for MRC, and `.spi` for SPIDER. Use separate working directories or unique
output names to preserve results from repeated runs.

## Limitations

This is a simplified Fourier-projection simulator, not a multislice
electron-scattering model for heavy-atom specimens. It generates boxed
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
