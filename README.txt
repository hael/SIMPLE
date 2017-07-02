
       _______ _____ _______  _____         _______
       |______   |   |  |  | |_____] |      |______
       ______| __|__ |  |  | |       |_____ |______  v2.5
 
simplecryoem.com

ABOUT

Single-particle IMage Processing Linux Engine is a program package for 
cryo-EM image processing, focusing on ab initio 3D reconstruction of 
single-particles with any point-group symmetry. The SIMPLE back-end 
consists of an object-oriented numerical library written in modern 
Fortran. The SIMPLE front-end consists of many standalone, interoperable 
components developed according to the "Unix toolkit philosophy".

SIMPLE is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the 
Free Software Foundation, either version 3 of the license, or (at your 
option) any later version. SIMPLE is distributed with the hope that it 
will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
the GNU General Public License for more details.

INSTALLATION INSTRUCTIONS

System Requirements:

    - Linux: (we use Ubuntu 15.04 and above)
    - MacOSX: 10.10 and above
    - CMake 3.2 and above
    - FFTW 3.3 and above
    - GNU toolchain (gcc & gfortran) 4.9 to 5.4

Installation:

1. Create the directory in which you are going to install SIMPLE 
(referred to here as <simple_path>):

    $ mkdir <simple_path>

2. Unzip the SIMPLE 2.5 tarball in this directory (assuming you have 
downloaded the tarball in the <downloads> directory):
    
    $ mv <downloads>/SIMPLE2.5.tgz <simple path>
    $ cd <simple path>
    $ tar -xzf SIMPLE2.5.tgz

3. Create a directory for the build:

    $ cd simple2.5
    $ mkdir build
    $ cd build

4. Compile and install SIMPLE 2.5:

    $ cmake ../
    $ make -j install

This will install SIMPLE in the 'build' directory. If you wish to 
provide an alternative installation directory, substitute step 4
with:

    $ cmake -DCMAKE_INSTALL_PREFIX=<alternative directory> ../
    $ make -j install

Step 4 assumes that gcc/gfortran and FFTW are installed in fairly 
standard directories on your machine. In case you have a more exotic 
setup you can provide the paths pointing to your custom gcc/gfortran & 
FFTW by substituting step 4 with:

    $ FC=<gcc/gfortran path> FFTW_DIR=<FFTW path> cmake ../
    $ make -j install

For instance, on MacOS 
- Macports users may use: FC=/opt/local/bin/gfortran FFTW_DIR=/opt/local;
- Fink users: FC=/sw/bin/gfortran FFTW_DIR=/sw/; and
- Homebrew users: FC=/usr/local/bin/gfortran FFTW_DIR=/usr/local/

5. To run SIMPLE, the bin and scripts paths need to be in the PATH 
environment variable and the lib path in the LD_LIBRARY_PATH variable. 
The SIMPLE_PATH environment variable must also be defined. The shell 
scripts add2.bashrc and add2.tcshrc containing the necessary instructions 
were generated during the build step. For immediate use for running and 
testing, execute

    $ source add2.bashrc

or, for TCSH/CSH users:

    $ source add2.tcshrc

For permanent installation BASH users should add the contents of add2.bashrc 
to your <HOME>/.bashrc

    $ cat add2.bashrc >> ~/.bashrc

or for TCSH/CSH users:

    $ cat add2.tcshrc >> ~/.tcshrc

To test the build, please execute

    $ make test
    $ ctest --output-on-failure

TESTING THE MOST IMPORTANT FEATURES

To ensure that SIMPLE has been correctly installed, we recommend 
running the application simple_test_install. It will test the most 
important components in the SIMPLE library  (those used by prime2D
and  prime3D. Execute

    $ simple_test_install 

The program will create its own folder SIMPLE_TEST_INSTALL*date*
where temporary files and information about each test are stored. Upon 
succesful completion you should see

    $ **** SIMPLE_TEST_INSTALL NORMAL STOP ****

simple_test_install can be executed anywhere. After execution, the folder 
created can be safely removed. If any of the individual tests fail an 
error message will be displayed. If you detect an error, please carefully
check the SIMPLE and FFTW installations and the gfortran version. If you 
still have issues, please file a help ticket on the webpage.
