#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script create_simulData.csh generates simulated data       !"
echo "!                for a given data set here polimerase                   !"
echo "!               Uses Prime2 and simulator from SIMPLE lib               !"
echo "!  How to use: csh create_simulData.csh ./polii.spi ./ 8 15 8 120 1.77  !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."

set model_in = $1            #Model in for consideration
set rundir_in = $2           #local dir
set nspace_in = $3           #number of simulated particles
set lp_in = $4               #low pass limit
set nthr_in = $5             #number of theads for the OpemMP 
set boxSize_in = $6          #box size
set smpd_in = $7             #smapling distance

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

set model = $model_in
set rundir = $rundir_in
set nspace = $nspace_in
set lp = $lp_in
set nthr = $nthr_in
set boxSize = $boxSize_in
set smpd = $smpd_in

set dir_name = $nspace$undscr"lp"$lp$undscr"bz"$boxSize$undscr"smpd"$smpd"_Bunch"
              printf "Creating data dir: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; mkdir $rundir$dir_name
              printf "Copying the model: "; tput bold;
tput setaf 1; printf "$model "; tput sgr0;
              printf "to directory: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; cp $model $rundir$dir_name

cd $dir_name
pwd

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Generating the simulated data...                           -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-            Generating random orientations simple_orisops...           -"
tput sgr0;
simple_orisops nptcls=$nspace outfile=rndoris.txt trs=3 ctf=yes kv=300 cs=2 defocus=2 deferr=1
tput setaf 2;
echo "-            Projecting the volume with simple_projvol...               -"
tput sgr0;
simple_projvol vol1=$model smpd=$smpd outstk=projs.spi oritab=rndoris.txt ctf=no
tput setaf 2;
echo "-            Projecting the volume with simple_simemimgs...             -"
tput sgr0;
simple_simemimgs vol1=$model nptcls=$nspace smpd=$smpd nthr=$nthr msk=50 snr=0.1 ctf=yes kv=300 cs=2 oritab=rndoris.txt outfile=simdoc.txt outstk=simimgs.spi stk=projs.spi bfac=50
tput setaf 2;
echo "-            Stack operation with simple_stackops...                    -"
tput sgr0;
simple_stackops smpd=$smpd stk=simimgs.spi outstk=simimgs_flip.spi oritab=simdoc.txt ctf=flip kv=300 cs=2
tput setaf 2;
echo "-            Reconstruct the volume with simple_eo_recvol...            -"
tput sgr0;
simple_eo_recvol smpd=$smpd nthr=$nthr stk=simimgs_flip.spi msk=50 oritab=simdoc.txt ctf=yes kv=300 cs=2 oritab=simdoc.txt
tput setaf 2;
echo "-            Stack operation with simple_stackops...                    -"
tput sgr0;
simple_stackops smpd=$smpd lp=$lp stk=simimgs_flip.spi shalgn=yes trs=3
tput setaf 2;
echo "-            Recenter the box with simple_resrange...                   -"
tput sgr0;
simple_resrange smpd=$smpd box=$boxSize
tput setaf 2;
echo "-            simple_rndrec...                                           -"
tput sgr0;
simple_rndrec nspace=$nspace smpd=$smpd nthr=$nthr stk=simimgs.spi msk=50
tput setaf 2;
echo "-            Prime 2...                                                 -"
tput sgr0;
simple_prime2 smpd=$smpd lp=$lp nthr=$nthr stk=simimgs_flip.spi vol1=recvol_state1.spi msk=50 oritab=rndoris.txt
tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

