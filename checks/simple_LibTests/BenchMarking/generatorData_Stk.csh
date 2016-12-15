#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7` # $8` #$9 $10 $11 $12 $13`
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
echo "!        1.77 #smpd value                                               !"
echo "!        8    #nspace                                                   !"
echo "!        8    #nthr                                                     !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set up_1dir = "../"
set undscr = "_"
set dot = "."
set ext_asc = "asc"
set ext_mrc = "mrc"
set spc = " "

set rundir_in     = $1  #local dir
set pdb_in        = $2  #input pdb
set mask_in       = $3  #maximum iterations
set boxSize_in    = $4  #box size
set smpd_in       = $5   #smapling distance
set nspace_in     = $6   #number of simulated particles
set nthr_in       = $7   #number of theads for the OpemMP 

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Maximum iterations: "; tput bold;
tput setaf 1; printf "$mask_in\n"; tput sgr0;
              printf "Input pdb particle: "; tput bold;
tput setaf 3; printf "$pdb_in\n"; tput sgr0;
              printf "BoxSize pixel: "; tput bold;
tput setaf 6; printf "$boxSize_in\n"; tput sgr0;
              printf "Sampling distance: "; tput bold;
tput setaf 8; printf "$smpd_in\n"; tput sgr0;
              printf "N particles to be simulated: "; tput bold;
tput setaf 7; printf "$nspace_in\n"; tput sgr0;
              printf "Number of threads: "; tput bold;
tput setaf 5; printf "$nthr_in \n"; tput sgr0;

set rundir     = $rundir_in
set pdb        = $pdb_in
set mask       = $mask_in
set boxSize    = $boxSize_in
set smpd       = $smpd_in
set nspace     = $nspace_in
set nthr       = $nthr_in

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Setting up the environment variables...                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

set mol_str = '1WCM'
set bz_str = 'bz'
set msk_str = 'msk'
set smpd_str = 'smpd'
set outstk_base = 'npart1_nptcls'
set outfile_base = 'oritab.txt'

set dir_base = $msk_str$mask$undscr$mol_str$boxSize
pwd
set out_dir = "$pdb$sptr$dir_base"
echo $out_dir
mkdir $out_dir

#"1WCM_bz128_smpd1.77.mrc"
set vol1 = "$pdb_in$sptr$mol_str$undscr$bz_str$boxSize$undscr$smpd_str$smpd$dot$ext_mrc"
echo $vol1

set bzsmpdmsk = "$bz_str$boxSize$undscr$smpd_str$smpd$undscr$msk_str" #'bz128_smpd1.77_msk'
set end_ext = "$mol_str$boxSize$dot$ext_mrc" #'1WCM128.mrc'

set outstk = "$out_dir$sptr$outstk_base$nspace$undscr$bzsmpdmsk$mask$undscr$end_ext"
set outfile = "$out_dir$sptr$outfile_base"
echo $outstk
echo $outfile

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Generating the different stackss...                        -"
echo "                   $dir_base                                             "
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

set prg="simimgs"
simple_exec prg=$prg vol1=$vol1 outstk=$outstk outfile=$outfile smpd=$smpd msk=54 nptcls=$nspace snr=.5 ctf=no nthr=$nthr ndiscrete=1000

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
