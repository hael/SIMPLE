       _______ _____ _______  _____         _______
       |______   |   |  |  | |_____] |      |______
       ______| __|__ |  |  | |       |_____ |______  v3b1

https://github.com/hael/SIMPLE/releases/tag/v3.0.0

PUBLIC UNIT TESTS
https://zenodo.org/records/18663655

ABOUT

SIMPLE is an end-to-end image processing and reconstruction platform. It takes raw electron microscopy data 
(movies, micrographs, particle images), processes and analyzes them through many stages (motion correction, 
CTF estimation, picking, classification, alignment, reconstruction), and produces scientifically meaningful 
outputs such as 2D class averages, 3D reconstructions, symmetry analysis, and atomic-level models. SIMPLE 
supports both interactive batch workflows and highly automated streaming pipelines, and it is designed to 
scale from a single workstation to distributed, high-performance computing environments. SIMPLE is a modular, 
high-performance scientific framework for cryo-EM and nanoparticle image analysis that combines low-level 
numerical kernels, rich domain abstractions, and high-level workflow orchestration into a scalable, test-
driven application suite. It is not a library in the narrow sense, nor a single executable—it is a full 
scientific computing environment, engineered for correctness, performance, extensibility, and real-world 
experimental workflows.

SIMPLE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the license, or (at your option) any later version. SIMPLE is distributed with the hope that it will be useful, but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

INSTALLATION INSTRUCTIONS

System Requirements:

    - Linux: (we use Ubuntu 16.04 and above)
    - MacOSX: 10.10 and above
    - CMake 3.5 and above
    - FFTW 3.3 and above
    - GNU toolchain (gcc & gfortran) 14.2 and above
    - libTIFF 4+ (for TIFF format support, ensure dev packages installed)
    - jbigkit 2+ (for TIFF format support, ensure dev packages installed)

Installation:

1. Obtain the source code for the desired SIMPLE release. 

    https://github.com/hael/SIMPLE/releases/tag/v3.0.0

    $ gunzip SIMPLE-3.0.0.tar.gz
    $ tar -xvf SIMPLE-3.0.0.tar
    $ cd SIMPLE-3.0.0

    Or to use the nightly build:

    $ git clone https://github.com/hael/SIMPLE.git
    $ cd SIMPLE

2. create a directory for the build

    $ mkdir build
    $ cd build

3. Compile and install in the 'build' directory

    $ cmake ..
    $ make -j install

For users who do not wish to install TIFF support (only required for motion correction), type instead:

    $ cmake -DUSE_LIBTIFF=OFF ..
    $ make -j install

Upon completion, compilation will finish with the following message

Installation complete.
==========================================================================
Please ensure the following variables are set properly in add2.*rc file:
    SIMPLE_EMAIL SIMPLE_QSYS SIMPLE_PATH
To use SIMPLE, append the relevant add2.* to your HOME shell rc file:
  bash$ cat add2.bashrc >> ~/.bashrc
  tcsh$ cat add2.tcshrc >> ~/.tcshrc
==========================================================================
Alternatively, for minimal installation to work correctly add:
<simple_path>/build/bin and
<simple_path>/build/scripts
to your PATH environment variable.
==========================================================================

4. To complete installation the relevant environment variables need to be set. To run SIMPLE the <simple_path>/build/bin and <simple_path>/build/scripts paths need to be in the PATH environment variable. The SIMPLE_PATH environment variable must also be defined. Shell scripts are provided for this purpose. BASH users should append the contents of add2.bashrc to their home ~/.bashrc:

    $ cat add2.bashrc >> ~/.bashrc

and for TCSH/CSH users:

    $ cat add2.tcshrc >> ~/.tcshrc

INSTALLATION NOTES

If you wish to provide an alternative installation directory, substitute step 3 with

    $ cmake -DCMAKE_INSTALL_PREFIX=<alternative directory> ..
    $ make -j install

Step 3 assumes that GNU gfortran and FFTW are installed in fairly standard directories on your machine. In case you have a more exotic setup you can provide the paths pointing to your custom gcc/gfortran & FFTW by substituting step 3 with

    $ FC=<gfortran path> CC=<gcc path> FFTW_DIR=<FFTW path> cmake ..
    $ make -j install

For instance, on MacOS such locations could be:

 - Macports:   FC=/opt/local/bin/gfortran FFTW_DIR=/opt/local
 - Fink users: FC=/sw/bin/gfortran FFTW_DIR=/sw/
 - Homebrew:   FC=/usr/local/bin/gfortran FFTW_DIR=/usr/local/

For more advanced FFTW installations, see *Installation of FFTW* below.

TESTING THE BUILD

To ensure that SIMPLE has been correctly installed, we recommend running the application simple_test_install. It will perform elementary tests of the base components in the SIMPLE library. Execute the following in a separate terminal to ensure the environment variables have been correctly set:

    $ simple_test_install

The program will create its own folder SIMPLE_TEST_INSTALL*date* where temporary files are stored. Upon succesful completion you should see

    $ **** SIMPLE_TEST_INSTALL NORMAL STOP ****

simple_test_install can be executed anywhere and the folder created can be safely removed. If any of the individual tests fail an error message will be displayed. If you detect an error, please carefully check the SIMPLE and FFTW installations and the gfortran version. If you still have issues, please file a help ticket on the webpage.

### Installation of FFTW

SIMPLE requires three versions of the FFTW library. A double-precision, a single-precision and a threaded-single precision build.

Copy the source from the FFTW homepage (www.fftw.org) and unpack it:
$ wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz
$ tar zxvf fftw-3.3.8.tar.gz
$ cd fftw-3.3.8

By default, the bootstrap script configures and builds the double precision and threaded libraries:
$ ./bootstrap.sh --prefix=/usr/local
$ make -j
$ sudo make install
$ make distclean

Now build the single-precision libraries:

$ ./bootstrap.sh --prefix=/usr/local --enable-single
$ make -j
$ sudo make install
$ make distclean

************************************************************************
INSTALLING ON A MAC THAT HASN'T HAD ANY SOFTWARE INSTALLED ON IT BEFORE


	1.	Install homebrew
$/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
Add brew to your path
	2.	Install fftw and libraries
$curl -o fftw-3.3.8.tar.gz http://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz
$tar zxvf fftw-3.3.8.tar.gz
build double precision and threaded libraries:
$cd fftw-3.8.8
$./bootstrap.sh –prefix=/usr/local 
$make -j
$sudo make install
$make distclean
Now build-single precision libraries:
$./bootstrap.sh –prefix=/usr/local –enable-single
$make -j
$sudo make install
$make distclean
$cd ..
	3.	Install git, gcc, libtiff, jbigkit and cmake – make sure you have a python version higher than 3.10 
$brew install git
$brew install cmake
$brew install gcc
$brew install libtiff
$brew install jbigkit
$brew install python@3.10
	4.	Pull and install SIMPLE
$git clone https://github.com/hael/SIMPLE.git
$cd SIMPLE
$mkdir build
$cd build
$cmake -D NICE=yes ..
$make -j install
$cat add2.bashrc >> ~/.zshrc
	5.	Start the SERVER – username and password are in the log to window
$source add2.bashrc
$nice_local

--------------------------------
TO UPDATE THE CODE

	1.	Go to the SIMPLE directory and clear out local changes
$cd ~/SIMPLE
$git stash
$git stash clear
	2.	Pull new code
$git pull
	3.	Install new code
$cd build
$make -j install
	4.	Start the SERVER again – username and password are in the log to window
$source add2.bashrc
$nice_local

*********************************************************************************



