#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4` # $5 $6 $7 $8 $9 $10 $11 $12 $13`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script parse_simulData.csh generates simulated data        !"
echo "!                for a given data set                                   !"
echo "!               Uses Prime3D and simulator from SIMPLE lib              !"
echo "! Usage: csh generate_data.csh {Input}                                  !"
echo "! Input: ./ /scratch/frederic/1WCM/1WCM.pdb  {...}                      !"
echo "! {...}: 8  #maxits                                                     !"
echo "!        1WCM #molecule                                                 !"
echo "!        70   #mask value                                               !"
echo "!        004  #lp value                                                 !"
echo "!        1.77 #smpd value                                               !"
echo "!        8    #nspace                                                   !"
echo "!        8    #nthr                                                     !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."
set ext_asc = "asc"
set ext_mrc = "mrc"
set spc = " "

set rundir_in     = $1  #local dir
set pdb_in        = $2  #input pdb
set mask_in       = $3  #maximum iterations
set boxSize_in    = $4  #box size

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Maximum iterations: "; tput bold;
tput setaf 1; printf "$mask_in\n"; tput sgr0;
              printf "Input pdb particle: "; tput bold;
tput setaf 3; printf "$pdb_in\n"; tput sgr0;

set rundir     = $rundir_in
set pdb        = $pdb_in
set mask       = $mask_in
set boxSize    = $boxSize_in


#set molecules  = $molecules_in
#set nspace     = $nspace_in
#set lp         = $lp_in
#set nthr       = $nthr_in
#set boxSize    = $boxSize_in
#set smpd       = $smpd_in


tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Setting up the environment variables...                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;


tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Reference volume generation...                             -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
#e2pdb2mrc.py --box 128,128,128 --apix 1.77  $pdb 1WCM_bz128_smpd1.77.mrc   #mask: 50 55
#e2pdb2mrc.py --box 160,160,160 --apix 1.416 $pdb 1WCM_bz160_smpd1.416.mrc  #mask: 65 70
#e2pdb2mrc.py --box 256,256,256 --apix 1.416 $pdb 1WCM_bz256_smpd1.416.mrc  #mask: 65 70 75 80 85 90 95 100 110 115
#e2pdb2mrc.py --box 288,288,288 --apix 1.011 $pdb 1WCM_bz288_smpd1.011.mrc  #mask: 90 95 100 105 110 115 120 125 130 135

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Generating the different volumes...                        -"
echo "-                       1WCM_bz128_msk                                  -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

















tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
