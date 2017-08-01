#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7 $8`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script create_simulData_clutering.csh generates simulated  !"
echo "!            data for given data set here polimerase                    !"
echo "! Use: csh create_simulData_clutering.csh ./polii.spi ./ 8 15 1 120 1.77 50  !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."

set model_in   = $1  #Model in for consideration
set rundir_in  = $2  #local dir
set nspace_in  = $3  #number of simulated particles
set lp_in      = $4  #low pass limit
set nthr_in    = $5  #number of theads for the OpemMP
set boxSize_in = $6  #box size
set smpd_in    = $7  #sampling distance
set mask_in    = $8  #Mask

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "Input nptcls: "; tput bold;
tput setaf 1; printf "$nspace_in\n"; tput sgr0;
              printf "low pass filter: "; tput bold;
tput setaf 2; printf "$lp_in \n"; tput sgr0;
              printf "mask: "; tput bold;
tput setaf 3; printf "$mask_in\n"; tput sgr0;
              printf "Sampling distance: "; tput bold;
tput setaf 4; printf "$smpd_in\n"; tput sgr0;
              printf "Model in for consideration: "; tput bold;
tput setaf 5; printf "$model_in\n"; tput sgr0;
              printf "Number of threads (OpenMP): "; tput bold;
tput setaf 6; printf "$nthr_in \n"; tput sgr0;
              printf "BoxSize: "; tput bold;
tput setaf 8; printf "$boxSize_in\n"; tput sgr0;

set lp      = $lp_in
set mask    = $mask_in
set smpd    = $smpd_in
set nthr    = $nthr_in

printf "Current directory: "; tput bold;
tput setaf 2; pwd

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                   Generating the simulated data...                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-                   Generating/testing clustering...                    -"
tput sgr0;
tput setaf 2;
echo "-                        dummy orientations                             -"
echo "-                   generating random orientations...                   -"
tput sgr0;
set prg="orisops"
simple_exec prg=$prg nptcls=2 ndiscrete=10 outfile=rndoris.txt trs=2
cp rndoris.txt tmp1.txt
cat rndoris.txt tmp1.txt tmp1.txt tmp1.txt tmp1.txt tmp1.txt tmp1.txt > tmp2.txt
rm tmp1.txt
mv tmp2.txt rndoris.txt
tput setaf 2;
echo "-                          reprojections...                             -"
tput sgr0;
set prg="projvol"
simple_exec prg=$prg vol1=$model_in smpd=$smpd outstk=projs.spi oritab=rndoris.txt nthr=$nthr
tput setaf 2;
echo "-                        initial clusters...                            -"
tput sgr0;
set prg = "prime2D_init"
simple_exec prg=$prg stk=projs.spi smpd=$smpd msk=$mask ncls=2 nthr=$nthr lp=20
tput setaf 2;
echo "-                        Actual clustering...                           -"
tput sgr0;
set prg = "prime2D"
simple_exec prg=$prg stk=projs.spi refs=startcavgs.spi smpd=$smpd_in oritab=prime2D_startdoc.txt ncls=2 msk=$mask_in lp=20 trs=2
tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0
