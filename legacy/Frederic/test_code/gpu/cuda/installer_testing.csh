#!/bin/bash
set ARGV = `basename -a $1 $2`
set sptr = "/"
set undscr = "_"
set dot = "."


tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!  Script install_MacOSX_binarties.csh installs SIMPLE librray binaries !"
echo "!  and dependencies:                                                    !"
echo "!                     1. gcc-4.9-bin.tar.gz                             !"
echo "!                     2. gfortran-4.9-bin.tar.gz                        !"
echo "!                     3. fftw-3.3.4.tar.gz                              !"
echo "!  Usage: csh install_MacOSX_binarties.csh {target path ex: ./ ./Simple}!"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput setaf 4;
echo "!                        WARNING...                                     !"
echo "!                                                                       !"
echo "!   This script will install GNU chain tool, gcc, g++ version 4.9:      !"
echo "!   cpp  g++  gcc  gcc-ar  gcc-nm  gcc-ranlib  gcov  gfortran           !"
echo "!             in: /usr/local/bin/                                       !"
echo "!   It wil also configure, compile and install fftw-3.4.4 for           !"
echo "! float, double, long-double and quad-precision with threads and OpenMp !"
echo "!             in: /usr/local/lib/                                       !"
echo "!                                                                       !"
echo "!   This will replace previous version previously installed             !"
echo "       Do you wish to proceed: yes or no?"; set input = $<
switch ("$input")
case *"yes"
    tput setaf 2;
    echo "you have entered: $input"
    echo "proceeding with installation..."
    breaksw
case *"no"
    tput setaf 1;
    echo "you have entered: $input, stop."
    exit
    breaksw
endsw
tput sgr0;

#Previous stuff
tput setaf 2;
echo "-                      fftw-3.3.4 library ...                           -"
tput sgr0;
tput setaf 6;
ls /usr/lib/x86_64-linux-gnu/*fftw*
tput sgr0;

tput setaf 2;
echo "-    GNU compilers /usr/local/bin/{gcc-4.9,g++-4.9, gfortran-4.9 ...    -"
tput sgr0;
tput setaf 5;
gcc-4.9 --version
tput setaf 2;
g++-4.9 --version
tput setaf 3;
gfortran-4.9 --version
tput sgr0;

set arbse = ("$1" "$2")

unset list
foreach dir (./)
foreach ext ("c" "f90" "tar")
 find $dir"/" -maxdepth 1 -name "*.$ext"
end
end

set extract_dir = ./latemp/latralala
set simple_home = "{full path to: }"$extract_dir$sptr"MacOSX_binaries"
echo "export SIMPLE_HOME=$simple_home                                          "

set i = 0
while ($i < $#arbse)
echo $i
echo $arbse[$i]

 @ i = $i + 1
end


foreach type ("float") # "double" "long-double" "quad-precision")

#type threads and OpenMP
tput setaf 2;
echo "-            Configuring the fftw-3.3.4 library                         -" 
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
#./configure CC=gcc-4.9 --enable-$type --enable-threads --enable-openmp
echo  "--enable-$type --enable-threads --enable-openmp"
tput setaf 2;
echo "-            compiling the fftw-3.3.4 library for                       -"
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
#make
tput setaf 2;
echo "-            installing the fftw-3.3.4 library                          -"
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
#make install
tput setaf 2;
echo "-    fftw-3.3.4 library for $type threads and OpenMP:Installed          -"
tput sgr0;
tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                       Checking installation ...                       -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
#ls /usr/lib/x86_64-linux-gnu/*fftw*



end
