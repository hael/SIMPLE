#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script parse_simulData.csh generates simulated data        !"
echo "!                for a given data set                                   !"
echo "!               Uses Prime3D and simulator from SIMPLE lib              !"
echo "! Usage: csh parse_simulData.csh {Input}                                !"
echo "! Input: ./ trpv1 0008 004 8 256 1.2156 70 yes no no 0 10               !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."
set ext = "asc"
set spc = " "
set file_comb = "z_prime3D"

set rundir_in     = $1   #local dir
set molecules_in  = $2   #File_Stack in for consideration
set nspace_in     = $3   #number of simulated particles
set lp_in         = $4   #low pass limit
set nthr_in       = $5   #number of theads for the OpemMP 
set boxSize_in    = $6   #box size
set smpd_in       = $7   #smapling distance
set mask_in       = $8   #Mask
set use_gpu_in    = $9   #use gpu yes or no
set bench_gpu_in  = $10  #bench gpu yes or no
set fix_gpu_in    = $11  #bench gpu yes or no
set set_gpu_in    = $12  #bench gpu yes or no
set maxits_in     = $13  #maximum iterations

tput bold;
tput setaf 4; printf "Input coming into create_simulData.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
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
tput setaf 1; printf "$maxits_in\n"; tput sgr0;

set rundir     = $rundir_in
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

              printf "Echoing data dir: "; tput bold;
tput setaf 2; printf "$rundir$dir_name\n"; tput sgr0;

cd $dir_name
printf "Files in current directory: "; tput bold;
tput setaf 2; pwd ;tput sgr0;
#listing the files in the current directory
ls
#collecting and parsing the information from the benchmarking files
#and putting the information in a common file
set strg_loop = ""
echo $maxits
set n=101
set n_max = 0
@ n_max = $n + $maxits 
while ( $n < $n_max )
      @ n = $n + 1
      echo $n
      set strg_loop = "$strg_loop$spc$n"
end
echo $strg_loop

set filename_base = "z_prime3D_iter"
#foreach iter (1002 1003 1004 1005 1006 1007 1008 1009 1010 1011)
foreach iter ($strg_loop)
set filename = $filename_base$iter$dot$ext

#cat $filename

set strg = ""

@ count = 0 #counter init
foreach LINE (`cat $filename`)
@ count ++
#echo $LINE
#echo $count
set s = $LINE
set strg = "$strg$spc$LINE"
end #end of the LINE loop

echo "gpu"$spc$use_gpu$spc"bench"$spc$bench_gpu$spc"fix"$spc$fix_gpu$spc"set"$spc$set_gpu$spc$strg$spc$molecules >> ../"gpu"$undscr$use_gpu$undscr"bench"$undscr$bench_gpu$undscr"fix"$undscr$fix_gpu$undscr"set"$undscr$set_gpu$undscr$file_comb$dot$ext

end #end of the iter loop 
