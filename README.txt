
       _______ _____ _______  _____         _______
       |______   |   |  |  | |_____] |      |______
       ______| __|__ |  |  | |       |_____ |______  v3_dev
 
simplecryoem.com

ABOUT

Single-particle IMage Processing Linux Engine is a program package for cryo-EM
image processing, focusing on ab initio 3D reconstruction of single-particles
with any point-group symmetry. The SIMPLE back-end consists of an
object-oriented numerical library written in modern Fortran. The SIMPLE
front-end consists of many standalone, interoperable components developed
according to the "Unix toolkit philosophy".

SIMPLE is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the license, or (at your option) any later
version. SIMPLE is distributed with the hope that it will be useful, but WITHOUT
ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

INSTALLATION INSTRUCTIONS

System Requirements:

    - Linux: (we use Ubuntu 15.04 and above)
    - MacOSX: 10.10 and above
    - CMake 3.2 and above
    - FFTW 3.3 and above
    - GNU toolchain (gcc & gfortran) 4.9 to 5.4

Installation:

1. Create the directory in which you are going to install SIMPLE (referred to
here as <simple_path>):

    $ mkdir <simple_path>

2. Unzip the SIMPLE 2.5 tarball in this directory (assuming you have downloaded
the tarball in the <downloads> directory):
    
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

This will install SIMPLE in the 'build' directory, clean out all unnecessary
files and will finish with the following message (a reminder for step 5 below):
"""
Installation complete.
==========================================================================
Please ensure the following variables are set properly in add2.*rc file:
    SIMPLE_EMAIL SIMPLE_QSYS SIMPLE_PATH 
To use SIMPLE, append the relevant add2.* to your HOME shell rc file:
  bash$ cat add2.bashrc >> ~/.bashrc
  tcsh$ cat add2.tcshrc >> ~/.tcshrc
==========================================================================
For minimal installation to work correctly add:
<your src path>/Simple-release/build/bin and
<your src path>/Simple-release/build/scripts
to your PATH environment variable.
==========================================================================
"""

When the build and installation directories are the same (default) and you are
happy with the install, you may want to clean compilation-generated and
unnecessary build files using distclean.

    $ make distclean

If you wish to provide an alternative installation directory, substitute step 4 with:

    $ cmake -DCMAKE_INSTALL_PREFIX=<alternative directory> ../
    $ make -j install

Step 4 assumes that GNU gfortran and FFTW are installed in fairly standard
directories on your machine. In case you have a more exotic setup you can
provide the paths pointing to your custom gcc/gfortran & FFTW by substituting
step 4 with:

    $ FC=<gfortran absolute path> FFTW_DIR=<FFTW path> cmake ../
    $ make -j install

For instance, on MacOS
 - Macports users may use: FC=/opt/local/bin/gfortran FFTW_DIR=/opt/local
 - Fink users: FC=/sw/bin/gfortran FFTW_DIR=/sw/ and
 - Homebrew users" FC=/usr/local/bin/gfortran FFTW_DIR=/usr/local/

5. Set the environment variables

To run SIMPLE the bin and scripts paths need to be in the PATH environment
variable. The SIMPLE_PATH environment variable must also be defined. The example 
shell scripts add2.bashrc and add2.tcshrc with the necessary instructions were 
generated during the build step.

For immediate use for running and testing:
    $ source add2.bashrc
or, for TCSH/CSH users:
    $ source add2.tcshrc

For permanent installation BASH users should add the contents of add2.bashrc to
your <HOME>/.bashrc:
    $ cat add2.bashrc >> ~/.bashrc
and for TCSH/CSH users:
    $ cat add2.tcshrc >> ~/.tcshrc

TESTING THE BUILD

To ensure that SIMPLE has been correctly installed, we recommend running the
application simple_test_install. It will test the most important components in
the SIMPLE library (those used by prime2D and prime3D). Execute the following in
a separate terminal to ensure the environment variables are set by your rc file:

    $ simple_test_install 

The program will create its own folder SIMPLE_TEST_INSTALL*date*
where temporary files and information about each test are stored. Upon 
succesful completion you should see

    $ **** SIMPLE_TEST_INSTALL NORMAL STOP ****

simple_test_install can be executed anywhere. After execution, the folder
created can be safely removed. If any of the individual tests fail an error
message will be displayed. If you detect an error, please carefully check the
SIMPLE and FFTW installations and the gfortran version. If you still have
issues, please file a help ticket on the webpage.


### Using the GUI in multiuser mode

The Simple GUI is designed to have a single instance of the server running on an accessible node, to which multiple users can connect using individual or shared login credentials.

In this mode, the Simple GUI server will run as the user who started the server. Thus, this user must have read and write access to the directories where processing is to be undertaken.

In order to use this mode, an htpasswd file must first be generated. The easiest way to do this is using the "htdigest" command from the apache2 package. Use the following command to create users:
	
	htdigest $SIMPLE_PATH/www/.htpasswd simple <username>

You will be prompted to enter the new user's password twice.

The Simple GUI server can now be started with the -m flag for multiuser mode. We recommend that this process is run in the background using the nohup command:

	nohup simple -m &



