       _______ _____ _______  _____         _______
       |______   |   |  |  | |_____] |      |______
       ______| __|__ |  |  | |       |_____ |______  v3b1

simplecryoem.com
https://github.com/hael/SIMPLE/releases/tag/v3.0.0

ABOUT

Single-particle IMage Processing Linux Engine (SIMPLE) is an open-source software package for analysis of cryogenic transmission electron microscopy (cryo-EM) movies of single-particles (Single-Particle Analysis, SPA). SIMPLE 3.0 provides many new features, including a full SPA processing pipeline and a streaming SPA platform for informing data acquisition in real time, using only minimal CPU computing resources. Our stream SPA tool implements the steps of anisotropic motion correction, particle identification and 2D clustering with automatic class rejection. SIMPLE 3.0 additionally features an easy-to-use web-based graphical user interface that can be run on any device (workstation, laptop, tablet or iphone) and supports remote multi-user execution over the network. The new project-based execution model automatically records the executed workflow, facilitates meta-data handling and greatly simplifies usage. Using SIMPLE 3.0, it is possible to automatically obtain a clean SP data set amenable to high-resolution 3D reconstruction directly upon completion of the data acquisition, without the need for extensive batch image processing.

SIMPLE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the license, or (at your option) any later version. SIMPLE is distributed with the hope that it will be useful, but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

INSTALLATION INSTRUCTIONS

System Requirements:

    - Linux: (we use Ubuntu 16.04 and above)
    - MacOSX: 10.10 and above
    - CMake 3.2 and above
    - FFTW 3.3 and above
    - GNU toolchain (gcc & gfortran) 5.5
    - libTIFF 4+ (for TIFF format support, ensure dev packages installed)
    - jbigkit 2+ (for TIFF format support, ensure dev packages installed)

Installation:

1. Obtain the source code for the desired SIMPLE release. 

    https://github.com/hael/SIMPLE/releases/tag/v3.0.0

    $ gunzip SIMPLE-3.0.0.tar.gz
    $ tar -xvf SIMPLE-3.0.0.tar
    $ cd SIMPLE-3.0.0

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

### Graphical User Interface (GUI)

Two versions of the GUI are installed during SIMPLE compilation: single-user and multi-user. Both are identical in appearance and functionality, but differ in the way in which they are accessed.

SINGLE-USER GUI

The single-user GUI is launched with the command 'simple'. The resulting window is visible only to the user on the current display. All data associated with the gui is stored in the users' directory, and jobs run as that user. When launched like this, the GUI is not accessible remotely. This GUI mode is recommended for evaluating the software, when running on a laptop or workstation (without remote access required), or in situations when users may not be able to remotely access the GUI e.g. due to network restrictions.

MULTI-USER GUI

The multi-user GUI is designed for a single instance of the server to be running on an accessible node, to which multiple users can connect using individual or shared login credentials. 

In this mode, the SIMPLE GUI server will run as the user who started the server. Thus, this user must have read and write access to the directories where processing is to be undertaken, as well as sufficient submission rights to any available cpu resources in a cluster. 

In order to use this mode, a file must be created at $SIMPLE_PATH/gui_data/users.htpasswd containing per-user login information in the form <USERNAME>:<PASSWORD>. It is recommended to change the permissions on this file so that only the user which runs the GUI server can read this file e.g. chmod 600 $SIMPLE_PATH/gui_data/users.htpasswd. As many users can be added as desired; each user will only see their projects. 

All GUI associated data is stored in $SIMPLE_PATH/gui_data/ so this directory may want to be backed up regularly. 

Once the login information has been entered, the GUI server can be started using the command 'simple_multiuser'. Alternatively, the process can be run in the background using the command 'nohup simple_multiuser &'. 

The GUI can now be accessed from any modern web browser by navigating to http://<address of node on which server is running>:8095. 

To facilitate users connecting externally to the cluster, it is possible to place the server behind a HTTP/HTTPS(recommended) proxy server such as apache2 or nginx. More information can be found at https://httpd.apache.org/docs/2.4/howto/reverse_proxy.html and https://docs.nginx.com/nginx/admin-guide/web-server/reverse-proxy/.

In the event that local firewall or security restictions on the cluster do not allow access via a HTTP/HTTPS proxy server, ssh may be used to forward traffic to the GUI server. The syntax of the ssh command is:
	ssh -L 8095:<name/address of machine running the SIMPLE GUI in the cluster>:8095 <your username>@<cluster login address>

For example:
	ssh -L 8095:simple_gui.madeupcluster.com:8095 madeupuser@madeupclusterlogin.com

The user can then use a browser on the machine that they SSD'd from to navigate to localhost:8095.
 

### Tutorials

The extensive SIMPLE manual describing the functionality of every back-end component has been replaced by the largely self-explanatory GUI that comes with embedded tutorials. When you launch the GUI you see three clickable icons in the bottom left corner. Clickling the rightmost of these - the lamp - opens the SIMPLE tutorials. The data set for the tutorials is available for download at http://simplecryoem.com/SIMPLE3.0/download.html

### Installation of FFTW

Simple requires three versions of the FFTW library. A double-precision, a single-precision and a threaded-single precision build.

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
