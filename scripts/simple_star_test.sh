!#/usr/bin/bash

SIMPLE_PATH=~/el85_scratch/eagerm/SIMPLE3.0/build
module load git cmake fftw/3.3.5-gcc5 relion/2.1 fftw/3.3.5-gcc5 gctf/1.06_cuda8 chimera eman/2.2
export PATH=$SIMPLE_PATH/bin:$PATH
cd ~/el85_scratch/stars_from_matt

mkdir simple_convert_star
cd simple_convert_star
## make new project

simple_exec prg=new_project projname=micrographs_all_gctf
cd micrographs_all_gctf
simple_exec prg=importstar_project starfile=../../micrographs_all_gctf.star

cd ..
simple_exec prg=new_project projname=Micrographs
cd Micrographs
simple_exec prg=importstar_project starfile=../../Micrographs/

## cryosparc
simple_exec prg=new_project projname=cryosparc_movies
cd cryosparc_movies
echo "../../cryosparc_exp000325_002.mrc" > filetab.txt
simple_exec prg=import_movies filetab=filetab.txt cs=2.7 kv=300 fraca=0.1 smpd=1


