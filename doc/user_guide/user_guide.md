# User Guide

`simple_exec` is the key SIMPLE executable, which has various 'commanders' to accomplish different tasks and run workflows available via the `prg=` keyword. The list of all commanders in the SIMPLE can be inspected with:
```shell
simple_exec prg=list
```

## Workflow Guides

- [reproject](reproject.md) — generate projections from a volume.
- [simulate_movie](simulate_movie.md) — generate a movie and its aligned
  average micrograph from a projection stack, with a PDB-to-movie tutorial.
- [simulate_particles](simulate_particles.md) — generate boxed particle
  images and their simulation metadata from a volume.
- [Simulated validation workflows](simulated_workflows.md) — run the
  built-in molecular workflows and quantitative checks with `simple_test_exec`.

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
simple_exec prg=print_project_info projfile=my_proj.simple
```
`projfile` is required; this command does not automatically select a project.
To inspect a project in a numbered execution directory, provide its path:

```shell
simple_exec prg=print_project_info projfile=1_import_particles/my_proj.simple
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


## Select a micrograph subset from an oriented project

Select exactly 100 micrographs using existing 3D particle scores and angular
coverage:

```shell
simple_exec prg=select_micrographs projfile=my_project.simple nmics=100 pgrp=c1
```

The command uses active particles in structural `state=1` by default. Change
`state=` to select another state. The default `subset_score=corr` uses the
last polar alignment score; `subset_score=corr_cart` selects the Cartesian
score instead. These scores are separate ranking proxies, not calibrated
probabilities. `subset_score=none` gives each eligible particle unit support.
Coverage selection requires updated particles with finite orientations and
complete existing project half-set provenance. It never creates a new split.

The default output is `my_project_mic_subset.simple` beside the input. Use
`outfile=another_subset.simple` for a different destination. The input and
image stacks are not modified. Output keeps original row indices and retained
particle poses, structural labels and halves; excluded rows are deselected.
Inherited output products, state weights and registered sigma state are
cleared, so reconstruct the selected particles before assessing map quality.

The output project has accompanying `.micrographs.csv` and `.coverage.csv`
reports containing selected micrographs, eligible particle counts, direction
coordinates and support before/after selection. `nspace=64` controls the
symmetry-aware direction grid and must be even. Set the correct molecular
symmetry with `pgrp=`.

The optimizer balances quality-weighted particle yield and directional support
using a logarithmic objective. It starts with a deterministic greedy subset,
then searches globally. `subset_max_nodes=100000` limits search work. The log
and micrograph report explicitly say whether the optimum of this objective
was proved; otherwise the subset is the best solution found within that
budget. This score does not guarantee an improved reconstructed volume.
This operation uses `mkdir=no`; choose an output path instead of an execution
directory. Existing output projects or reports are never overwritten.

To maximize accepted particle count without using orientations or scores:

```shell
simple_exec prg=select_micrographs projfile=my_project.simple nmics=100 subset_mode=particles
```

## Select a direction-balanced particle subset

Select exactly 100000 refined particles using the best score within each
symmetry-aware viewing-direction group:

```shell
simple_exec prg=select_particles projfile=my_project.simple nsample=100000 nspace=64 pgrp=c1
```

Use the reconstruction's symmetry in `pgrp`. `nspace` must be positive and
even. The default `state=1` restricts selection to that structural state.
The default `subset_score=corr` ranks the last polar alignment scores;
use `subset_score=corr_cart` for Cartesian refinement or polishing scores.
Higher scores rank first. Use comparable final-pass scores; this command
does not rescore particles, and fractional updates can leave older scores.

Quotas are equal across occupied direction bins until small bins are
exhausted, then redistributed. The final partial quota round retains the
highest available scores, giving exactly `nsample` particles. The command
rejects requests exceeding the eligible population. Inactive particles,
other states, unrefined rows and nonfinite poses or scores are excluded.
Complete existing half-set provenance is required. Half labels are rebalanced
by default; use `eo_rebalance=no` to preserve the existing labels.

The command requires `mkdir=yes` and creates a numbered execution directory
containing the modified copy of the input project, `selected_particles.mrc`,
and `selected_particles.txt`. The project copy is pruned: project particle
row N, stack image N and text row N describe the same retained particle.
Logical stack rows retain their microscope metadata but all reference the
single generated stack. Original `pind`, poses, shifts and states are
preserved. The source project and input image stacks are unchanged.

The default `eo_rebalance=yes` creates a new particle-level split balanced
by direction and score. Within every viewing bin, adjacent
score-ranked particles are paired into opposite halves, and unpaired
particles balance the total counts. Each bin and the total population
differ by at most one particle between halves. Pair orientation greedily
balances score sums within the bin. Both project particle segments and
the exported table receive the new labels; stack provenance records a
per-particle split so later processing preserves them.

Rebalancing splits micrograph groups and mixes the previously
refined halves. It does not preserve the original gold-standard
independence. Do not reuse the old half-maps or FSC as validation of the
new partition, and be aware that correlated neighbouring picks can now
fall in opposite halves. The source project's half tags are unchanged.

The text file uses native SIMPLE `key=value` orientation records, including
Euler angles, shifts, defocus, astigmatism, phase shift, voltage, Cs,
amplitude contrast and sampling. Units remain native: shifts in pixels,
defocus in micrometres, Euler/astigmatism angles in degrees, phase in radians,
voltage in kV, Cs in mm, and sampling in Angstroms per pixel.
Inherited map outputs and registered sigma state are cleared from the
subset project. Reconstruct the subset before judging map quality; score
selection alone does not guarantee improved resolution.

## NICE (New Interface for CryoEM)
Launch the Graphical User Interface (GUI) with:
```shell
nice_local
```

Start the SERVER – username and password are in the log to window

Coming Soon ...
