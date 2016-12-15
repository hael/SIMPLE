#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7 $8 $9 $10`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script create_simulData.csh generates simulated data       !"
echo "!                for a given data set here polimerase                   !"
echo "!               Uses Prime3D and simulator from SIMPLE lib              !"
echo "! Usage: csh create_simulData.csh {Input}                               !"
echo "! Input: ./polii.spi ./ 8 5 8 240 1.77 50 yes/no yes/no                 !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."

set model_in     = $1   #Model in for consideration
set rundir_in    = $2   #local dir
set nspace_in    = $3   #number of simulated particles
set lp_in        = $4   #low pass limit
set nthr_in      = $5   #number of theads for the OpemMP 
set boxSize_in   = $6   #box size
set smpd_in      = $7   #smapling distance
set mask_in      = $8   #Mask
set use_gpu_in   = $9   #use gpu yes or no
set bench_gpu_in = $10  #bench gpu yes or no

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Input model: "; tput bold;
tput setaf 1; printf "$model_in\n"; tput sgr0;
              printf "N particles to be simulated: "; tput bold;
tput setaf 3; printf "$nspace_in\n"; tput sgr0;
              printf "low pass filter: "; tput bold;
tput setaf 4; printf "$lp_in \n"; tput sgr0;
              printf "Number of threads: "; tput bold;
tput setaf 5; printf "$nthr_in \n"; tput sgr0;
              printf "BoxSize pixel: "; tput bold;
tput setaf 6; printf "$boxSize_in\n"; tput sgr0;
              printf "Sampling distance: "; tput bold;
tput setaf 8; printf "$smpd_in\n"; tput sgr0;
              printf "mask: "; tput bold;
tput setaf 4; printf "$mask_in\n"; tput sgr0;
              printf "Use GPU ---> yes/no: "; tput bold;
#tput setaf 8; printf "$use_gpu_in\n"; tput sgr0;
switch ("$use_gpu_in")
case *"yes"
    tput setaf 2;
    printf "you have entered: $use_gpu_in ---> "
    printf "proceeding with Hybrid GPU CPU calculation...\n"; tput sgr0;
    breaksw
case *"no"
    tput setaf 1;
    printf "you have entered: $use_gpu_in ---> "
    printf "proceeding with CPU calculation only...\n"; tput sgr0;
    breaksw
endsw
printf "Use GPU ---> yes/no: "; tput bold;
#tput setaf 8; printf "$use_gpu_in\n"; tput sgr0;
switch ("$bench_gpu_in")
case *"yes"
    tput setaf 2;
    printf "you have entered: $bench_gpu_in ---> "
    printf "proceeding with Hybrid GPU CPU calculation...\n"; tput sgr0;
    breaksw
case *"no"
    tput setaf 1;
    printf "you have entered: $bench_gpu_in ---> "
    printf "proceeding with CPU calculation only...\n"; tput sgr0;
    breaksw
endsw

set model = $model_in
set rundir = $rundir_in
set nspace = $nspace_in
set lp = $lp_in
set nthr = $nthr_in
set boxSize = $boxSize_in
set smpd = $smpd_in
set mask    = $mask_in
set use_gpu = $use_gpu_in
set bench_gpu = $bench_gpu_in

set dir_name = "gpu"$undscr$use_gpu$undscr"bench"$undscr$bench_gpu$undscr$nspace$undscr"lp"$lp$undscr"bz"$boxSize$undscr"smpd"$smpd$undscr"msk"$mask"_Bunch"

              printf "Creating data dir: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; mkdir $rundir$dir_name
              printf "Copying the model: "; tput bold;
tput setaf 1; printf "$model "; tput sgr0;
              printf "to directory: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; cp $model $rundir$dir_name

cd $dir_name
printf "Current directory: "; tput bold;
tput setaf 2; pwd

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Generating the simulated data...                           -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-            Generating random orientations simple_orisops...           -"
tput sgr0;
set prg="makeoris"
simple_exec prg=$prg nptcls=$nspace outfile=rndoris.txt trs=0 ctf=yes kv=300 cs=2 defocus=2 deferr=1
tput setaf 2;
echo "-            Projecting the volume with simple_projvol...               -"
tput sgr0;
set prg="projvol"
simple_exec prg=$prg vol1=$model smpd=$smpd outstk=projs.spi oritab=rndoris.txt ctf=no nthr=$nthr
tput setaf 2;
echo "-            Projecting the volume with simple_simemimgs...             -"
tput sgr0;
#ctf=yes kv=300 cs=2
set prg="simimgs"
simple_exec prg=$prg vol1=$model nptcls=$nspace smpd=$smpd nthr=$nthr msk=$mask snr=0.5 oritab=rndoris.txt outfile=simdoc.txt outstk=simimgs.spi
tput setaf 2;
echo "-            Stack operation with simple_stackops...                    -"
tput sgr0;
set prg="ctfops"
simple_exec prg=$prg smpd=$smpd stk=simimgs.spi outstk=simimgs_flip.spi oritab=simdoc.txt ctf=flip kv=300 cs=2
tput setaf 2;
echo "-            Reconstruct the volume with simple_eo_recvol...            -"
tput sgr0;
set prg="eo_recvol"
simple_exec prg=$prg smpd=$smpd nthr=$nthr stk=simimgs_flip.spi msk=$mask oritab=simdoc.txt ctf=flip kv=300 cs=2 oritab=simdoc.txt
tput setaf 2;
echo "-            Stack operation with simple_stackops...                    -"
tput sgr0;
set prg="stackops"
simple_exec prg=$prg smpd=$smpd lp=$lp stk=simimgs_flip.spi shalgn=yes trs=3
tput setaf 2;
echo "-            Recenter the box with simple_resrange...                   -"
tput sgr0;
set prg="resrange"
simple_exec prg=$prg smpd=$smpd box=$boxSize
tput setaf 2;
echo "-            Prime3D init...                                            -"
tput sgr0;
set prg="prime3D_init"
simple_exec prg=$prg nspace=$nspace smpd=$smpd nthr=$nthr stk=simimgs_flip.spi msk=$mask lp=$lp
tput setaf 2;
echo "-            Prime3D...                                                 -"
tput sgr0;
set prg="prime3D"
simple_exec prg=$prg smpd=$smpd lp=$lp nthr=$nthr stk=simimgs_flip.spi vol1=startvol_state01.spi msk=$mask oritab=prime3D_startdoc.txt use_gpu=$use_gpu bench_gpu=$bench_gpu nspace=$nspace
tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

