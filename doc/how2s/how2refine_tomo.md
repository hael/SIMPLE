S=~/src/SIMPLE/scripts/aet_mat2simple.py; D=~/Downloads/TomoFromPeter
mkdir -p ~/fept_tomo && cd ~/fept_tomo
python3 $S stack  $D/FePt_Denoised_Projections.mat Mproj projs.mrcs --smpd 0.3725 --pad 320
python3 $S angles $D/FePt_raw_projections_and_angles.mat Angles tilts.txt
simple_exec prg=new_project projname=fept && cd fept
single_exec prg=import_tilt_series projfile=fept.simple stk=../projs.mrcs angtab=../tilts.txt angconv=tilt smpd=0.3725 mkdir=no
single_exec prg=refine_tomo projfile=fept.simple mskdiam=108 nthr=10