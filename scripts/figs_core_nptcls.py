# * Generation of figures and movies of the core nanoparticles
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#
# Figures must be generated using UCSF Chimera on one of the designated Linux desktop systems where the Chimera graphical interface is available.
# Prior to execution, update the following variables in this script:
#
#     # Trajectory file paths (directory containing trajectories data)
#     path      = ["/home/meanapanedar2/elmlund/NanoX/new_data/docked_map_analysis"]
#     # Trajectory naming convention (name to be used in the manuscript)
#     traj_name = ["new_growth"]
#     # Time-segment identifier (naming convention for trajectory subdirectories)
#     segment   = "stage"
#     # Smpd Pixel size in Angstrom
#     smpd      = 0.264
#     # atomic element
#     element   = "Pt"
#
# Script Descriptions
# ───────────────────
#
#     gen_core_figs.py - Generates all figures of the core of the nanoparticles.
#
# Execution Commands
# ───────────────────
#
# Run the scripts using the following commands:
#
# python  gen_core_figs.py
#!/usr/bin/env python3

from PIL import Image, ImageOps, ImageDraw, ImageFont
import os
import re
from pathlib import Path

# Trajectory file paths (directory containing trajectories data)
path      = ["/home/meanapanedar2/elmlund/NanoX/new_data/docked_map_analysis",
             "/home/meanapanedar2/elmlund/NanoX/new_data/docked_map_analysis"]
# Trajectory naming convention (name to be used in the manuscript)
traj_name = ["trajectory_1",
             "trajectory_2"]
# Time-segment identifier (naming convention for trajectory subdirectories)
segment   = "stage"
# Smpd Pixel size in Angstrom
smpd      = 0.264
# atomic element
element   = "Pt"

# recreate trajs dir
if os.path.exists("trajs"):
    os.system("rm -rf trajs")
    os.system("rm -rf *.mp4")
os.mkdir("trajs")
os.chdir("trajs")
total_heigh  = 2000 * len(path)
core_img     = [ "_ATMS_core_b1_fringe_b0.pdb", "ATMS_core.pdb" ] 
core_img1    = [ "_ATMS_fringe.pdb", "ATMS_core.pdb" ] 
# Run SINGLE core atoms analysis 
for i in range(len(path)):
    print(">>> RUNNING SINGLE_CORE_ATOMS_ANALYSIS for ",traj_name[i])
    # make / enter per-trajectory folder
    os.mkdir(traj_name[i])
    os.chdir(traj_name[i])
    # write list of pdb files into pdbfiles.txt
    cmd_ls = f'ls -v "{path[i]}/{segment}"*/*ATMS.pdb > pdbfiles.txt'
    os.system(cmd_ls)
    # run SINGLE (uses pdbfiles.txt)
    cmd_single = f'single_exec prg=core_atoms_analysis pdbfiles=pdbfiles.txt smpd={smpd} element={element}'
    os.system(cmd_single)
    # List ATMS core files
    os.chdir("1_core_atoms_analysis")
    # generate png files
    print(">>> GENERATING IMAGES ")
    with open("chimera_script.py", "w") as f:
        f.write("from chimera import runCommand as rc\n")
        f.write("import sys\n")
        f.write("rc(\"windowsize 2000 2000\")\n")
        f.write("rc(\"windoworigin 0 0\")\n")
        cnt=0 
        for c in range(len(core_img)):
            # core_b1_fringe_b0
            p = Path.cwd()
            sfile_list = sorted(
                f.name for f in p.iterdir()
                if f.is_file() and segment and core_img[c] in f.name
            )
            png_list = [str(Path(f).with_suffix(".png")) for f in sfile_list]
            for file in reversed(sfile_list):
                figure = str(Path(file).with_suffix(".png"))
                f.write('file = "{}"\n'.format(file))
                f.write('figure = "{}"\n'.format(figure)) 
                f.write("rc(\"open \" + file)\n")
                if cnt==0 :
                    f.write("rc(\"center\")\n")
                    f.write("rc(\"focus\")\n")
                f.write("rc(\"scale 1.0\")\n")
                f.write("rc(\"colordef tpink 1. 5. 7. .0\")\n")
                f.write("rc(\"setattr p color tpink \")\n")
                f.write("rc(\"represent sphere\")\n")
                f.write("rc(\"set projection ortho\")\n")
                f.write("rc(\"back solid white\")\n")
                f.write("rc(\"preset apply pub 3\")\n")
                f.write("rc(\"unset depthCue\")\n")
                f.write("rc(\"turn y  5 \")\n")
                f.write("rc(\"turn x  2\")\n")
                f.write("rc(\"rangecol bfactor,a key 0. #f0f921 1. #0d0887 \")\n")
                f.write("rc(\"copy file \" + figure + \" png \")\n")
                f.write("rc(\"close all\")\n")
                cnt+=1
        f.write("rc(\"stop now\")\n")
    os.system("chimera chimera_script.py")

    # generate movies 
    home_dir = os.path.expanduser("~")
    tmp_dir  = os.path.join(home_dir, "tmp")
    print(">>> GENERATING MOVIES")
    with open("chimera_movie_script.py", "w") as f:
        f.write("from chimera import runCommand as rc\n")
        f.write("import os\n")
        f.write("import re\n")
        f.write("os.system(\"rm -rf {}\")\n".format(tmp_dir))
        f.write("os.system(\"mkdir -p {}\")\n".format(tmp_dir))
        f.write("rc(\"windowsize 2000 2000\")\n")
        f.write("rc(\"windoworigin 0 0\")\n")
        cnt=0 
        for c in range(len(core_img1)):
            f.write('print("Making movie ----", "{}", "{}")\n'.format(traj_name[i], core_img1[c]))
            coret = core_img1[c].replace(".pdb", "")
            movie = traj_name[i] + "_" + coret + ".mp4"
            f.write('movie = "{}"\n'.format(movie)) 
            f.write("rc(\"close all\")\n")
            f.write("rc(\"reset\")\n")
            sfile_list = sorted(
                f.name for f in p.iterdir()
                if f.is_file() and segment and core_img1[c] in f.name
            )
            if cnt==0 : f.write("rc(\"open core.pdb\")\n")
            for s in range(len(sfile_list)):
                f.write("rc(\"open {}\")\n".format(sfile_list[s]))
            f.write("rc(\"focus\")\n")
            f.write("rc(\"scale 0.8\")\n")
            f.write("rc(\"colordef tpink 1. 5. 7. .0\")\n")
            f.write("rc(\"setattr p color tpink\")\n")
            f.write("rc(\"represent sphere\")\n")
            f.write("rc(\"set projection ortho\")\n")
            f.write("rc(\"back solid white\")\n")
            f.write("rc(\"preset apply pub 3\")\n")
            f.write("rc(\"unset depthCue\")\n")
            f.write("rc(\"rangecol bfactor,a key 0. #f0f921 1. #0d0887\")\n")
            f.write("rc(\"movie record directory  {}\")\n".format(tmp_dir))
            f.write("rc(\"turn y  5 \")\n")
            f.write("rc(\"turn x  2 \")\n")
            f.write('rc("savepos start_pos")\n')
            for s in range(len(sfile_list)):
                f.write('rc("~modeldisplay")\n')
                if cnt==0 : 
                    f.write('rc("modeldisplay #0")\n')
                    f.write("rc(\"modeldisplay #{}\")\n".format(s+1))
                else:
                    f.write("rc(\"modeldisplay #{}\")\n".format(s))
                f.write("rc(\"wait \")\n")
                f.write("rc(\"wait 45  \")\n")
                f.write("rc(\"turn y 1 45  \")\n")
                f.write("rc(\"wait 45  \")\n")
                f.write("rc(\"turn y -1 45  \")\n")
                f.write("rc(\"wait 45  \")\n")
                f.write("rc(\"turn y -1 45  \")\n")
                f.write("rc(\"wait 45  \")\n")
                f.write("rc(\"turn y 1 45  \")\n")
                f.write("rc(\"wait 45  \")\n")
                f.write("rc(\"reset start_pos  \")\n")
                f.write("rc(\"wait 45  \")\n")
            f.write('rc("movie stop")\n')
            f.write("rc(\"movie encode output ../../../{} wait true quality highest\")\n".format(movie))
            f.write("os.system(\"rm -f {}/*\")\n".format(tmp_dir))
            cnt+=1
        f.write('rc("close all")\n')
        f.write('rc("close session")\n')
        f.write('rc("stop now")\n')
    os.system("chimera chimera_movie_script.py")

    # generate full figure
    for c in range(len(core_img)):
        # core_b1_fringe_b0
        p = Path.cwd()
        sfile_list = sorted(
            f.name for f in p.iterdir()
            if f.is_file() and segment and core_img[c] in f.name
        )
        png_list = [str(Path(f).with_suffix(".png")) for f in sfile_list]
        row          = 0
        total_width  = 2000 * len(sfile_list)
        figure       = Image.new('RGB', (total_width,2000), color=(255,255,255))
        col          = 0
        for im in png_list:
            im   = Image.open(im)
            figure.paste(im, (col,row))
            col+=2000
        draw = ImageDraw.Draw(figure)
        figfile = "../../" + traj_name[i] + "_" + core_img[c] + ".png"
        figure.save(figfile)

    # combine all figures 
    total_width  = 2000 * len(sfile_list)
    total_height = 4000 
    figure       = Image.new('RGB', (total_width,total_height), color=(255,255,255))
    col = 0
    row = 0
    for c in range(len(core_img)):
        im   = "../../" + traj_name[i] + "_" + core_img[c] + ".png"
        im   = Image.open(im)
        figure.paste(im, (col,row))
        row+=2000
    draw = ImageDraw.Draw(figure)
    figfile = "../../" + traj_name[i] + ".png"
    figure.save(figfile)
    
    # go back up for next iteration
    os.chdir("../..")



