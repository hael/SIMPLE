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

## mini_stream and abinitio2d are similar workflows the former is 

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

## Run Abinitio2D workflow
```
simple_exec prg=abinitio2D ncls=100 mskdiam=180 nthr=36
```

## Open the file 7_abinitio2D/cavgs_iter0*_ranked.mrc (last iteration) using e2display.py from EMAN2
cd 4_mini_stream
E2display.py

Click the middle mouse button to open up a popup in e2display
Select "Del" button in the popup window. Then select all the images that should be removed.

Next click the middle mouse button to bring up the popup again and click "Save" save the remaining images into a file named "selected.spi"

## Map class average selection
```
simple_exec prg=map_cavgs_selection stk2=./7_abinitio2D/selected.spi prune=yes
```

## Run Abinitio3D
```
simple_exec prg=abinitio3D pgrp=c1 mskdiam=180 nthr=40
```


# NICE (New Interface for CryoEM)
Launch the Graphical User Interface (GUI) with:
```shell
nice_local
```

Start the SERVER – username and password are in the log to window

Coming Soon ...