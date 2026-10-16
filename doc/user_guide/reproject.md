# reproject — Generate Projections from a Volume

`simple_exec prg=reproject` generates a stack of 2D projections from a 3D
density volume using Fourier interpolation. It can generate spiral projection
directions or reuse orientations from a SIMPLE orientation table.

## Inputs and setup

Use an installation configured as described in
[Installation](../installation.md), with `simple_exec` on `PATH`.
Supply a centered cubic molecular volume in MRC or SPIDER format. Set `smpd`
to its actual sampling in Å/pixel; this option does not resample the volume.

The examples assume a 192 × 192 × 192 MRC volume at 1.3 Å/pixel, with the
molecule fitting within a 180 Å mask. Replace `/absolute/path/reference.mrc`
with your volume and adjust the sampling, symmetry, and mask for your molecule.
The examples use `mkdir=no` inside fresh working directories.

## Examples

### Generate spiral directions

```shell
mkdir reproject_spiral
cd reproject_spiral

simple_exec prg=reproject \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 nspace=10 pgrp=c1 mskdiam=180 \
  outstk=projections.mrcs mkdir=no nthr=1
```

This writes ten 192 × 192 projection images and their orientation records.
Use the resulting clean stack as input to
[simulate_movie](simulate_movie.md) when preparing synthetic micrographs.

### Reuse supplied orientations

This example first creates a ten-record orientation table, saves it under a
different name, and reuses it for a second projection stack:

```shell
mkdir reproject_supplied
cd reproject_supplied

simple_exec prg=reproject \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 nspace=10 pgrp=c1 mskdiam=180 \
  outstk=initial_projections.mrcs mkdir=no nthr=1

cp reproject_oris.txt input_oris.txt

simple_exec prg=reproject \
  vol1=/absolute/path/reference.mrc \
  smpd=1.3 nspace=10 pgrp=c1 mskdiam=180 \
  oritab=input_oris.txt outstk=supplied_projections.mrcs mkdir=no nthr=1
```

For an existing SIMPLE text orientation table, replace `input_oris.txt` with
its path. For a single-volume run, use records assigned to `state=1`.
The current CLI requires `nspace` and `pgrp` even with `oritab`; once the
table is loaded, its record count determines the number of output images.
Keep the input table separate from the fixed `reproject_oris.txt` output.

## Controls and defaults

| Option | Meaning | Effective default or requirement |
| --- | --- | --- |
| `vol1` | Input density volume | Required |
| `smpd` | Sampling in Å/pixel | Supply the actual input sampling; standalone parameter initialization is 1.3 |
| `nspace` | Number of generated projection directions | Required by the CLI; with `oritab`, the table count is used |
| `pgrp` | Point group for generated directions | Required by the CLI; examples use `c1` |
| `oritab` | SIMPLE orientation table | Omitted: generate spiral directions |
| `mskdiam` | Soft-mask diameter in Å | Omitted: this commander does not mask the volume |
| `neg` | Invert projection contrast | `no` |
| `outstk` | Projection-stack filename | `reprojs.mrcs` |
| `mkdir` | Create a numbered execution directory | `no` |
| `nthr` | OpenMP thread count | Runtime OpenMP maximum when omitted; examples use 1 |

For generated directions, `pgrp` restricts sampling to the asymmetric unit.
Point-group families include cyclic (`cN`), dihedral (`dN`), tetrahedral
(`t`), octahedral (`o`), and icosahedral (`i`) symmetry. Supplied orientations
are read from the table rather than generated from this setting.

If `mskdiam` is supplied, SIMPLE converts it to a pixel radius and applies a
soft spherical mask before projection. A mask that is too large for the box
falls back to the box-derived mask with a warning. Choose the physical
diameter to include your molecule.

## Outputs

| File | Contents |
| --- | --- |
| `outstk` (default `reprojs.mrcs`) | Projection images, one per generated direction or input orientation record |
| `reproject_oris.txt` | SIMPLE orientation records corresponding to the projection stack |

The output filename extension selects the image format: `.mrc` or `.mrcs`
for MRC, and `.spi` for SPIDER. The default `reprojs.mrcs` is an MRC stack.
Output image dimensions follow the input volume box.

An existing `outstk` is deleted before the new stack is written.
`reproject_oris.txt` is a fixed output name and is rewritten on another run
in the same directory. Use separate working directories or preserve both
outputs before rerunning.

## Limitations

This command produces projections without simulated detector noise, movie
motion, or CTF effects. Use [simulate_movie](simulate_movie.md) for full-frame
movies or [simulate_particles](simulate_particles.md) for noisy boxed
particle images. Fourier projection is a simplified density-image model,
not a multislice electron-scattering simulation.

The examples cover a single molecular volume. Orientation-table state labels
matter: state-zero records produce blank images, and active records must
refer to a supplied volume. Changing `smpd` alone does not change the volume's
pixel dimensions or provide a physical resampling operation.

For built-in molecular fixtures and quantitative workflow checks, see
[Simulated validation workflows](simulated_workflows.md).
