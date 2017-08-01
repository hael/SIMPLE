#!/bin/bash
set ARGV = `basename -a $1` # $2`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            Script install_MacOSX_binarties.csh installs the binaries  !"
echo "!                and dependencies:                                      !"
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
tput setaf 4;
echo "!                                                                       !"
tput sgr0;

set sptr = "/"
set undscr = "_"
set dot = "."

#set rundir_in  = $1       #local dir
set extract_dir_in = $1  #extraction dirrectory dir

#set rundir_tmp =  $rundir_in
set extract_dir = $extract_dir_in

#cd $rundir_tmp
set rundir = `pwd`
printf "Current directory: "; tput bold;
tput setaf 2; pwd

tput bold;
tput setaf 4; printf "Input coming into install_MacOSX_binarties.csh:\n"; tput sgr0;
              printf "local exec dir: "; tput bold;
tput setaf 2; printf "$rundir\n"; tput sgr0;
              printf "Extraction directory: "; tput bold;
tput setaf 1; printf "$extract_dir_in\n"; tput sgr0;

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-            Installing the binaries for MacOSX...                      -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-     Creating the directory of extraction and untaring ball in it ...  -"
tput sgr0;
#MacOSX_binaries_Install.tar
mkdir -p $extract_dir && tar xf MacOSX_binaries_Install.tar -C $extract_dir
#gcc-4.9-bin.tar.gz
tput setaf 2;
echo "-            Installing gcc-4.9-bin.tar.gz compilers ...                -"
tput sgr0;
sudo tar xfz gcc-4.9-bin.tar.gz -C $extract_dir #TODO: replace with /usr/local/bin in sudo mode
#gfortran-4.9-bin.tar.gz
tput setaf 2;
echo "-            Installing gfortran-4.9-bin.tar.gz compilers ...           -"
tput sgr0;
sudo tar xfz gfortran-4.9-bin.tar.gz -C $extract_dir #TODO: replace with /usr/local/bin in sudo mode
tput setaf 2;
echo "-            Extracting fftw-3.3.4.tar.gz library ...                   -"
tput sgr0;
tar xfz fftw-3.3.4.tar.gz

#fftw-3.4.3
tput setaf 2;
echo "-            Moving to fftw-3.3.4.tar.gz directory ...                  -"
tput sgr0;
cd fftw-3.3.4
printf "Current directory: "; tput bold;
tput setaf 2; pwd
tput sgr0;

#Loop over types:
foreach type ("float" "single" "double" "long-double" "quad-precision")
#type threads and OpenMP
tput setaf 2;
echo "-            Configuring the fftw-3.3.4 library                         -" 
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
sudo ./configure CC=gcc-4.9 --enable-$type --enable-threads --enable-openmp
tput setaf 2;
echo "-            compiling the fftw-3.3.4 library for                       -"
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
sudo make
tput setaf 2;
echo "-            installing the fftw-3.3.4 library                          -"
echo "-              for $type threads and OpenMP ...                         -"
tput sgr0;
#sudo make install
tput setaf 2;
echo "-    fftw-3.3.4 library for $type threads and OpenMP:Installed          -"
tput sgr0;
tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                       Checking installation ...                       -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
ls /usr/lib/x86_64-linux-gnu/*fftw*
end #end foreach loop for type

#Doner and checkers.

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                              Done.                                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                       Checking installation ...                       -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;
tput setaf 2;
echo "-                                                                       -"
echo "-                      fftw-3.3.4 library ...                           -"
echo "-                                                                       -"
tput sgr0;
tput setaf 6;
ls /usr/lib/x86_64-linux-gnu/*fftw* #TODO: replace with /usr/local/lib/*
tput sgr0;

tput setaf 2;
echo "-                                                                       -"
echo "-    GNU compilers /usr/local/bin/{gcc-4.9,g++-4.9, gfortran-4.9} ...   -"
echo "-                                                                       -"
tput sgr0;
tput setaf 6;
#ls /usr/local/bin/g*  #TODO: uncomment the following
tput sgr0;
tput setaf 1;
gcc-4.9 --version      #TODO: replace with for MacOSX /usr/local/bin/gcc-4.9 --version
tput setaf 2;
g++-4.9 --version      #TODO: replace with for MacOSX /usr/local/bin/g++-4.9 --version
tput setaf 3;
gfortran-4.9 --version #TODO: replace with for MacOSX /usr/local/bin/gfortran-4.9 --version
tput sgr0;

tput bold; tput setaf 4;
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "-                    Finhing up the installation ...                    -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
tput sgr0;

cd $rundir
sudo rm -rf fftw-3.3.4
sudo chown -R $USER $extract_dir

printf "Finishing directory: "; tput bold;
tput setaf 2; pwd

tput bold; tput setaf 4;
set simple_home = "{full path i.e.: }"$extract_dir$sptr"MacOSX_binaries"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!            You now need to add the environment variable in            !"
echo "!                your ~/.bashrc or shell rc:                            !"
echo "!                                                                       !"
tput setaf 1;
echo "                             in bash                                     "
tput setaf 2;
echo "#Simple path (in bash)                                                   "
echo "export SIMPLE_HOME=$simple_home                                          "
echo 'export SIMPLEBIN=$SIMPLE_HOME/bin                                        '
echo 'export SIMPLEBINTESTS=$SIMPLE_HOME/bin/bin_tests                         '
echo 'export PATH=$PATH:$SIMPLEBIN:$SIMPLEBINTESTS:$SIMPLE_HOME/scripts        '
echo '                                                                         '
tput setaf 1;
echo "                     in shell or cshell or tshell                        "
tput setaf 2;
echo "#Simple path (in cshell)                                                 "
echo "setenv SIMPLE_HOME $simple_home                                          "
echo 'setenv SIMPLEBIN $SIMPLE_HOME/bin                                        '
echo 'setenv SIMPLEBINTESTS $SIMPLE_HOME/bin/bin_tests                         '
echo 'setenv PATH=(${SIMPLEBIN} ${SIMPLEBINTESTS} ${SIMPLE_HOME}/scripts $PATH '
tput setaf 1;
echo '                                                                         '
echo "                    source ~/.bashrc or source ~/.cshrc                  "
tput setaf 4;
echo "!                                                                       !"
echo "    Check ownership of installation path: ls -al $extract_dir           "
echo "!        if it has root and not {username} then change it by:           !"
echo "!                                                                       !"
tput setaf 2;
echo "                sudo chown -R {username} $extract_dir                    "
tput setaf 4;
echo "!                                                                       !"
echo "!          You may now start using SIMPLE from command line             !"
echo "      and go to: cd $extract_dir/MacOSX_binaries                         "
echo "!                       to launch the checks                             "
echo "!                                                                       !"
tput setaf 2;
echo "                        csh launch_checks.csh                            "
tput setaf 4;
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

