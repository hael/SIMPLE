# simulate_movie — Simulate a Movie and Micrograph

`simple_exec prg=simulate_movie` places projection images into a synthetic
micrograph and generates noisy movie frames with global drift and CTF
effects. It also writes an aligned average using the known simulated shifts.

## Inputs and setup

Use an installation configured as described in
[Installation](../installation.md), with `simple_exec` on `PATH`.
Supply a stack of square, clean particle projections. Each stack image is
placed once, so the stack length determines the number of particles.
[reproject](reproject.md) describes how to prepare this stack from a volume.
The simulation commands below do not require a SIMPLE project.

Set `smpd` to the input sampling in Å/pixel; it does not resample the images.
Use clean projections to avoid adding the noise and CTF from an earlier
particle simulation a second time.

## Simulating movies from a PDB structure

The workflow is **PDB coordinates → MRC volume → projection stack → movie**.
`simulate_movie` takes 2D projections rather than PDB coordinates or a volume.
If you already have a projection stack, start at
[step 3](#3-place-the-projections-into-a-simulated-movie) and adjust the sampling
and movie dimensions for your inputs.

The following example uses a molecule that fits within a 180 Å mask, a
128-pixel particle box and a sampling distance of 2 Å/pixel. Adjust the box
and mask to fit your structure. Run the commands in a new directory containing
your input `molecule.pdb`; `mkdir=no` keeps the outputs in that directory.

### 1. Convert Cartesian atomic coordinates to an MRC volume

```shell
simple_exec prg=pdb2mrc pdbfile=molecule.pdb \
  outvol=molecule.mrc pdbout=molecule_centered.pdb \
  smpd=2.0 vol_dim=128 center_pdb=yes mkdir=no
```

PDB Cartesian coordinates are in Ångströms. SIMPLE constructs an electrostatic
potential map from the atomic coordinates using Gaussian representations of
atomic scattering factors; it does not simply put one identical point at each
atom position.

| Parameter | Meaning |
| --- | --- |
| `smpd=2.0` | Sampling distance: 2 Å per voxel. |
| `vol_dim=128` | Cubic output dimensions: 128 × 128 × 128 voxels, spanning 256 Å along each axis in this example. |
| `center_pdb=yes` | Center the atomic model in the volume. |
| `outvol` | Output MRC volume. |
| `pdbout` | Output coordinates corresponding to the centered model. |

Choose an even box size with enough empty space around the molecule for all
projection directions. A smaller sampling distance gives a finer grid, but
requires more voxels to cover the same physical volume. If `vol_dim` is omitted,
SIMPLE determines the dimensions from the molecular extent.

### 2. Generate clean particle projections

```shell
simple_exec prg=reproject vol1=molecule.mrc outstk=projections.mrcs \
  smpd=2.0 nspace=100 pgrp=c1 mskdiam=180 neg=no nthr=1 mkdir=no
```

This creates 100 projections at reference directions distributed over the
orientation space. `pgrp=c1` assumes no symmetry; use the appropriate point
group if known. `mskdiam` is the soft-mask diameter in Ångströms and must enclose
the molecule without cutting it off. Keep `smpd` consistent with the volume.
See [reproject](reproject.md) for its controls and output orientation table.

### 3. Place the projections into a simulated movie

```shell
SIMPLE_SEED=20260923 simple_exec prg=simulate_movie stk=projections.mrcs \
  smpd=2.0 xdim=4096 ydim=4096 nframes=8 \
  snr=0.2 trs=2.0 kv=300 cs=2.7 fraca=0.1 \
  ctf=yes defocus=2.0 bfac=200 nthr=1 mkdir=no
```

The example creates a 4096 × 4096 movie with eight frames. The simulator
places the particles at random, non-overlapping positions; all frames share
this layout. See [Outputs](#outputs) for the movie, aligned average, and metadata.

A fixed `SIMPLE_SEED` and `nthr=1` help reproduce the random inputs with the
same SIMPLE build and environment; change the seed for a different placement,
motion and noise realization.

## Controls and defaults

These defaults describe the current implementation; some differ from values
displayed in command help.

| Option | Meaning | Effective default or requirement |
| --- | --- | --- |
| `stk` | Stack of square projection images to place | Supply an existing stack |
| `smpd` | Sampling in Å/pixel | Supply the actual input sampling; standalone parameter initialization is 1.3 |
| `xdim`, `ydim` | Output image dimensions in pixels | Set both explicitly; otherwise `xdim` is derived as half the particle box and `ydim` remains zero |
| `nframes` | Number of movie frames | 30 |
| `snr` | Per-frame signal-to-noise setting | Supply a positive value; initialized to zero |
| `trs` | Maximum random-walk increment per axis, in pixels | 3; cumulative drift can be larger |
| `kv` | Acceleration voltage in kV | 300 |
| `cs` | Spherical aberration in mm | 2.7 |
| `fraca` | Amplitude contrast fraction | 0.1 |
| `defocus` | Nominal underfocus in µm | 2 |
| `bfac` | CTF envelope B-factor in Å² | 200; fixed across frames rather than accumulated with dose |
| `mkdir` | Create a numbered execution directory | `yes`; the example uses `no` in a fresh directory |
| `nthr` | OpenMP thread count | Runtime OpenMP maximum when omitted; the example uses 1 |

For each movie, the two defocus components are randomly offset from
`defocus` by up to 0.2 µm, with a random astigmatism angle. One CTF is used
throughout that movie. The current generator applies CTF filtering even if
`ctf=no` is supplied; that option does not provide a CTF-free movie.

`trs` controls the random increments used to construct the drift trajectory,
not a bound on the total displacement. The actual shifts are recorded in
the output metadata.

## SNR and the movie noise model

SNR means **signal-to-noise variance ratio**, not an amplitude ratio or a value
in decibels. For a single Gaussian-noise addition, SIMPLE sets the noise
standard deviation to `sqrt(variance(input image) / snr)`. For example,
`snr=0.2` corresponds to five times as much noise variance as input-image
variance at that individual noise-addition step. Lower positive SNR values
give noisier images. Always supply a positive `snr` for movie simulation;
zero is not a noise-free setting.

The movie simulator uses two noise additions for each frame:

1. It shifts the contrast-inverted base image, normalizes it, and adds
   zero-mean Gaussian noise with an SNR of `snr / 0.2`.
2. It applies the CTF and B-factor envelope to both the signal and this first
   noise component. The filtered noise is spatially correlated; SIMPLE calls
   this component *pink noise*. It is not an explicitly generated `1/f` noise
   spectrum.
3. It normalizes the filtered image again and adds independent zero-mean white
   Gaussian detector noise with an SNR of `snr / 0.8`.

The 0.2 and 0.8 factors are nominal pre-CTF and detector-noise weights. Because
filtering and renormalization occur between the additions, they do not enforce
an exact 20%/80% split of the final noise variance. Likewise, the requested
`snr` is a per-frame simulation control, not an exact measured final-particle
SNR. The image variance includes empty background, so particle occupancy also
affects its interpretation. Averaging correctly aligned frames reduces the
independent noise; the average's SNR is therefore different from a frame's SNR.

These details describe `simulate_movie`; the separate
[simulate_particles](simulate_particles.md) command adds Gaussian detector
noise after optical filtering and uses a different noise model.

## Outputs

With MRC input as above, the simulator writes:

| File | Contents |
| --- | --- |
| `simulate_movie.mrc` | The stack of noisy, shifted movie frames |
| `optimal_movie_average.mrc` | One micrograph formed by undoing the known shifts and averaging |
| `simulate_movie_params.txt` | One record containing `dfx`, `dfy`, `angast`, `phshift`, and per-frame shifts `x1`, `y1`, etc. |

The aligned average still contains CTF effects and noise. Use it as a
synthetic integrated micrograph, or process `simulate_movie.mrc` through
motion correction first. The metadata records movie motion and CTF; it
does not export particle placement coordinates or per-particle orientations.
If the input was generated by `reproject`, its `reproject_oris.txt` describes
the projections rather than their positions in the micrograph.

The output basenames are fixed; `outstk` does not rename the movie output.
For multiple movies, use separate output directories or rename all three
outputs after each run. Reusing the same directory and filenames can
overwrite earlier results. Image output extensions follow the input format;
the filenames above are for MRC input.

## Placement and model limits

The occupancy check uses the full square particle boxes:

```text
occupancy = number_of_stack_images * box * box / (xdim * ydim)
```

The simulator stops above 0.55 occupancy. Each particle box must be smaller
than both micrograph dimensions. Random placement can also fail below the
limit if it cannot find enough non-overlapping positions. Reduce the stack
size (for example, reduce `nspace` when generating projections) or enlarge
the micrograph in that case. The example occupies approximately 9.8% of the
image area: `100 × 128² / 4096²`.

This is a simplified simulator for processing experiments and regression
checks. It uses a static projected scene, global translations, CTF filtering,
and Gaussian noise. It does not model Poisson electron counting, a calibrated
detector response, local motion, dose-dependent specimen damage, or multislice
electron scattering.

For built-in movie-to-3D and micrograph-to-2D checks, see
[Simulated validation workflows](simulated_workflows.md). For boxed particle
images generated directly from a volume, see
[simulate_particles](simulate_particles.md).
