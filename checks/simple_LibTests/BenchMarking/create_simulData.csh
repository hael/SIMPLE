#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15`
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

set rundir_in     = $1   #local dir
set file_stack_in = $2   #File_Stack in for consideration
set file_vol_in   = $3   #File_Stack in for consideration
set molecules_in  = $4   #File_Stack in for consideration
set nspace_in     = $5   #number of simulated particles
set lp_in         = $6   #low pass limit
set nthr_in       = $7   #number of theads for the OpemMP 
set boxSize_in    = $8   #box size
set smpd_in       = $9   #smapling distance
set mask_in       = $10  #Mask
set use_gpu_in    = $11  #use gpu yes or no
set bench_gpu_in  = $12  #bench gpu yes or no
set fix_gpu_in    = $13  #bench gpu yes or no
set set_gpu_in    = $14  #bench gpu yes or no
set maxits_in     = $15  #maximum iterations

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Input file_stack: "; tput bold;
tput setaf 1; printf "$file_stack_in\n"; tput sgr0;
              printf "Input file_vol: "; tput bold;
tput setaf 1; printf "$file_vol_in\n"; tput sgr0;
              printf "Input molecules: "; tput bold;
tput setaf 2; printf "$molecules_in\n"; tput sgr0;
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
printf "Bench GPU ---> yes/no: "; tput bold;
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
printf "Fix GPU ---> yes/no: "; tput bold;
#tput setaf 8; printf "$use_gpu_in\n"; tput sgr0;
switch ("$fix_gpu_in")
case *"yes"
    tput setaf 2;
    printf "you have entered: $fix_gpu_in ---> "
    printf "Fixing the GPU for the calculation...\n"; tput sgr0;
    breaksw
case *"no"
    tput setaf 1;
    printf "you have entered: $fix_gpu_in ---> "
    printf "proceeding without fixing the GPU...\n"; tput sgr0;
    breaksw
endsw
              printf "set GPU: "; tput bold;
tput setaf 1; printf "$set_gpu_in\n"; tput sgr0;
              printf "Maximum iterations: "; tput bold;
tput setaf 2; printf "$maxits_in\n"; tput sgr0;

set rundir     = $rundir_in
set file_stack = $file_stack_in
set file_vol   = $file_vol_in
set molecules  = $molecules_in
set nspace     = $nspace_in
set lp         = $lp_in
set nthr       = $nthr_in
set boxSize    = $boxSize_in
set smpd       = $smpd_in
set mask       = $mask_in
set use_gpu    = $use_gpu_in
set bench_gpu  = $bench_gpu_in
set fix_gpu    = $fix_gpu_in
set set_gpu    = $set_gpu_in
set maxits     = $maxits_in

set dir_name = "gpu"$undscr$use_gpu$undscr"bench"$undscr$bench_gpu$undscr"fix"$undscr$fix_gpu$undscr"set"$undscr$set_gpu$undscr$nspace$undscr"lp"$lp$undscr"bz"$boxSize$undscr"smpd"$smpd$undscr"msk"$mask$undscr$molecules

              printf "Creating data dir: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; mkdir $rundir$dir_name
#              printf "to directory: "; tput bold;
#tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0; cp $file_stack $rundir$dir_name

cd $dir_name
printf "Current directory: "; tput bold;
tput setaf 2; pwd

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Generating the simulated data...                           -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-            Prime3D...                                                 -"
tput sgr0;
set prg="prime3D"
simple_exec prg=$prg smpd=$smpd lp=$lp nthr=$nthr stk=$file_stack vol1=$file_vol msk=$mask use_gpu=$use_gpu bench_gpu=$bench_gpu fix_gpu=$fix_gpu set_gpu=$set_gpu nspace=$nspace maxits=$maxits
tput setaf 2;
echo "-            Gouping the results into a common file...                  -"
#TODO: need to fgix the maxit from above
csh ../parse_simulData.csh $rundir $molecules $nspace $lp $nthr $boxSize $smpd $mask $use_gpu $bench_gpu $fix_gpu $set_gpu $maxits

#removing unsed files for to avoid disk over filling
rm recvol_state01_iter*.mrc

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

