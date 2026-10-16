# User Guide

`simple_exec` is the key SIMPLE executable, which has various 'commanders' to accomplish different tasks and run workflows available via the `prg=` keyword. The list of all commanders in the SIMPLE can be inspected with:
```shell
simple_exec prg=list
```

## Creating a New Project
```
simple_exec prg=new_project projname=my_proj
```
This creates a my_proj directory and places a new SIMPLE binary file with the name my_proj.simple in the directory.

To create SIMPLE project file in the current directory use
```
simple_exec prg=new_project projname=my_proj dir=.
```

Similarly, to place the SIMPLE project file in a specific directory use the following:
```
mkdir my_dir
simple_exec prg=new_project projname=my_proj dir=path/to/directory
```
Note that the path to the dictory must exist, otherwise SIMPLE with raise an error that the directory does not exist.

## Print the Project Information
```
simple_exec prg=print_project_info
```
When run from within the project directory SIMPLE will automatically figure out the most recent project files to show.

If there are multiple project files in a directory, then provide the specific project file to print the info for:
```

```

## Create a list of absolute paths to EM movies in movies.txt from the directory ../movies/
filetab_movs.pl ../movies/

## Import  EM movies in .tiff format into the project file
```
simple_exec prg=import_movies cs=2.8 fraca=0.1 kv=200 smpd=0.885 filetab=movies.txt
```

## Motion Correction
```
simple_exec prg=motion_correct total_dose=52 nparts=2 nthr=20
```

## Estimate CTF
```
simple_exec prg=ctf_estimate nparts=2 nthr=40
```

## Create a list of absolute paths to .mrc files in filetab.txt from motion corrected micrographs (./2_motion_correct)
```
filetab_mrc.pl 2_motion_correct/
simple_exec prg=mini_stream cs=2.8 fraca=0.1 kv=200 smpd=0.885 nthr=36 filetab=filetab.txt
```

## mini_stream and solve2d are similar workflows the former is 

## Open the file 4_mini_stream/cavgs_iter0*_ranked.mrc (last iteration) using e2display.py from EMAN2
cd 4_mini_stream
E2display.py

Click the middle mouse button to open up a popup in e2display
Select "Del" button in the popup window. Then select all the images that should be removed. Next click "Save" top open the save window and save the remaining images into a file named "selected.spi"

## Convert the spi to mrc
```
simple_exec prg=convert stk=4_mini_stream/selected.spi outstk=selected.mrc smpd=0.885
```

## 2D Automasking
```
simple_exec prg=pick pickrefs=selected.mrc nparts=2 nthr=32
```

## Run Solve2D workflow
```
simple_exec prg=solve2D ncls=100 mskdiam=180 nthr=36
```

## Open the file 7_solve2D/cavgs_iter0*_ranked.mrc (last iteration) using e2display.py from EMAN2
cd 4_mini_stream
E2display.py

Click the middle mouse button to open up a popup in e2display
Select "Del" button in the popup window. Then select all the images that should be removed.

Next click the middle mouse button to bring up the popup again and click "Save" save the remaining images into a file named "selected.spi"

## Map class average selection
```
simple_exec prg=map_cavgs_selection stk2=./7_solve2D/selected.spi prune=yes
```

## Run Solve3D
```
simple_exec prg=solve3D pgrp=c1 mskdiam=180 nthr=40
```


## Simulating movies from a PDB structure

The workflow is **PDB coordinates → MRC volume → projection stack → movie**.
`simulate_movie` takes a stack of 2D particle projections, rather than a PDB
file or a 3D volume. These simulation commands do not require a SIMPLE project.

The following example uses a molecule that fits within a 180 Å mask, a
128-pixel particle box and a sampling distance of 2 Å/pixel. Adjust the box
and mask to fit your structure. Run the commands in a new directory containing
your input `molecule.pdb`; `mkdir=no` keeps the outputs in that directory.

### 1. Convert Cartesian atomic coordinates to an MRC volume

```bash
simple_exec prg=pdb2mrc pdbfile=molecule.pdb \
  outvol=molecule.mrc pdbout=molecule_centered.pdb \
  smpd=2.0 vol_dim=128 center_pdb=yes mkdir=no
```

PDB Cartesian coordinates are in Ångströms. SIMPLE constructs an electrostatic
potential map from the atomic coordinates using Gaussian representations of
atomic scattering factors; it does not simply put one identical point at each
atom position.

| Parameter | Meaning |
|---|---|
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

```bash
simple_exec prg=reproject vol1=molecule.mrc outstk=projections.mrcs \
  smpd=2.0 nspace=100 pgrp=c1 mskdiam=180 neg=no nthr=1 mkdir=no
```

This creates 100 projections at reference directions distributed over the
orientation space. `pgrp=c1` assumes no symmetry; use the appropriate point
group if known. `mskdiam` is the soft-mask diameter in Ångströms and must enclose
the molecule without cutting it off. Keep `smpd` consistent with the volume.
`reproject` supplies clean projections, leaving CTF filtering and noise addition
to the movie simulator.

### 3. Place the projections into a simulated movie

```bash
SIMPLE_SEED=20260923 simple_exec prg=simulate_movie stk=projections.mrcs \
  smpd=2.0 xdim=4096 ydim=4096 nframes=8 \
  snr=0.2 trs=2.0 kv=300 cs=2.7 fraca=0.1 \
  ctf=yes defocus=2.0 bfac=200 nthr=1 mkdir=no
```

Each image in the input stack is placed once at a random, non-overlapping
position. Frames share this particle layout and are related by global random
translations. `trs` controls the maximum random increment along each image
axis, in pixels; these increments accumulate over the frames.

The example creates a 4096 × 4096 movie with eight frames. `kv` is the accelerating
voltage in kV, `cs` the spherical aberration in mm, and `fraca` the amplitude
contrast fraction. `defocus` is the nominal underfocus in µm; the simulator
perturbs the two astigmatic defocus values around it and chooses an astigmatism
angle. `bfac` controls the CTF envelope damping in Å². These optical parameters
are shared across the movie's frames.

Keep particle occupancy low enough for placement to succeed. The simulator
checks the fraction `number of projections × box² / (xdim × ydim)` and rejects
values above 0.55; random placement may fail below that limit too. This example
occupies approximately 9.8% of the image area. Reduce `nspace` in the projection
step or increase the movie dimensions if placement fails.

Outputs are:

- `simulate_movie.mrc`: the movie frame stack.
- `optimal_movie_average.mrc`: an average aligned using the known simulated
  shifts. It still contains noise; it is not a noise-free reference.
- `simulate_movie_params.txt`: the simulated shifts and CTF parameters.

Use a separate directory for each movie. A fixed `SIMPLE_SEED` and `nthr=1`
help reproduce the random inputs with the same SIMPLE build and environment;
change the seed for a different placement, motion and noise realization.

### SNR and the movie noise model

SNR means **signal-to-noise variance ratio**, not an amplitude ratio or a value
in decibels. For a single Gaussian-noise addition, SIMPLE sets the noise
standard deviation to `sqrt(variance(input image) / snr)`. For example,
`snr=0.2` corresponds to five times as much noise variance as input-image
variance at that individual noise-addition step. Lower positive SNR values
give noisier images. Always supply a positive `snr` for movie simulation;
zero is not a noise-free setting.

The movie simulator uses two noise additions for each frame:

1. It shifts the clean base image, inverts particle contrast, normalizes the
   image, and adds zero-mean Gaussian noise with an SNR of `snr / 0.2`.
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

This is a simplified testing model: Gaussian noise and global translations,
without explicit Poisson electron counting, a calibrated detector response,
local motion or time-dependent radiation damage. The B-factor envelope is
fixed rather than accumulated with dose. These details describe
`simulate_movie`; the separate `simulate_particles` command currently adds
Gaussian detector noise after optical filtering and should not be assumed to
use the same two-component model.

# NICE (New Interface for CryoEM)
Launch the Graphical User Interface (GUI) with:
```shell
nice_local
```

Start the SERVER – username and password are in the log to window

Coming Soon ...
