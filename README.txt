Single-particle IMage Processing Linux Engine (SIMPLE) does ab initio
3D reconstruction, heterogeneity analysis, and high-resolution
refinement. The SIMPLE back-end consists of an object-oriented
numerical library with a single external dependency-the Fastest
Fourier Transform in the West (FFTW) (Frigo and Johnson, 2005). The
SIMPLE front-end consists of a few standalone, interoperable
components developed according to the `Unix toolkit
philosophy'. SIMPLE is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
license, or (at your option) any later version. This program is
distributed with the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details. You should have received a copy of the GNU General
Public License along with this program. If not, see
http://www.gnu.org/licenses/. 

----------------------------------
Updates and inclusion of new code:
----------------------------------

4th of March
-retructing of the library
27th of March
-Automation of the compilatiooj process
-Makefile_genAll.pl (uses simple_user_input.pl)

----------------------------------
System requirements
----------------------------------

Hardware

[CPU]
-MacOsX (yosemite and El Capitan)
-linux (fully supported Debians: Mint-17.1, Ubuntu (14.10, 15.04,
 15.10, 16.04 LTS), SUSE-13.2. Earlier versions are also supported but
 not tested)
-Intel Xeon5 2620 v2 and v3 but creates dependence on OS such as
 Ubuntu and hence the CUDA tool toolkit compability issues.  
[CPU/GPU]
-(Optional) a GPU card that supports CUDA (7.0 and 7.5) or OpenCL code on
  it. Here we have used NVIDIA K20c and K40.
-The supported achitecture of 3.5 and above but set as 3.5 to avoid
 conflicts with other cards

Software

[CPU]
- GNU chain tool 4.9 and above.
- Lapack and the Blas
- FFTW-3 (fast fourrier transforms)
[CPU/GPU]
- Only GNU chain tool 4.9
- CUDA development kit (7.0 and 7.5)
- NVIDIA_GPU_Computing_SDK/ (under Linux) GPU Computing (MacOsx)
  The two are the same different name.
- MAGMA library (magma-1.6.1)
- cuDNN v4 and v5
- OpenMPI preferably the MPI-3 but MPI-2 or above

----------------------------------
Instalation (MacOSX)
----------------------------------

first decide if you will use CPU [See (A)] or hybrid CPU/GPU
architecture [See (A) and (B)]




In the folder ./scripts there are some template for the input files
used to compile the SIMPLE, for example.

Template_MacOSX_10.9.5x64_CPU_simple_user_input.pm (also for 10.11.4)
Template_Mint_17.1x64_CPU_simple_user_input.pm
Template_SUSE_13.2x64_CPU_simple_user_input.pm
Template_FEDORA_21x64_CPU_simple_user_input.pm
Template_Ubuntu_14.10x64_CPU_simple_user_input.pm
Template_Ubuntu_15.10x64_CPU_simple_user_input.pm

You may use and Copy these for the one that suits your
platform into simple_user_input.pm file

>$ cp scripts/Template_MacOSX_10.9.5x64_CPU_simple_user_input.pm simple_user_input.pm

or you may edit simple_user_input.pm yourself using your favorite
editor.

------------
CPU
------------

(A) if you use only CPU, you only need to 

(1) change the path of the root of the application in simple_user_input.pm

$SIMPLE_PATH="{the full path of where simple is installed";

(2) If not already installed, install the FFTW-3 library this can be obtained from:

http://www.fftw.org/install/mac.html

Follow the instruction there but typically you will need to do a 

$./configure
$make
$sudo make install

and (you need to do it at least twice, for the floats double long and quad)

$./configure --enable-floats --enable-threads
$make
$sudo make install

check that you have in you lib folder (typically)

/usr/local/lib/libfftw3.a
/usr/local/lib/libfftw3.la
/usr/local/lib/libfftw3_threads.a
/usr/local/lib/libfftw3_threads.la
/usr/local/lib/libfftw3f.a
/usr/local/lib/libfftw3f.la
/usr/local/lib/libfftw3f_omp.a
/usr/local/lib/libfftw3f_omp.la
/usr/local/lib/libfftw3f_threads.a
/usr/local/lib/libfftw3f_threads.la
/usr/local/lib/libfftw3l.a
/usr/local/lib/libfftw3l.la
/usr/local/lib/libfftw3l_omp.a
/usr/local/lib/libfftw3l_omp.la
/usr/local/lib/libfftw3l_threads.a
/usr/local/lib/libfftw3l_threads.la
/usr/local/lib/libfftw3q.a
/usr/local/lib/libfftw3q.la
/usr/local/lib/libfftw3q_omp.a
/usr/local/lib/libfftw3q_omp.la
/usr/local/lib/libfftw3q_threads.a
/usr/local/lib/libfftw3q_threads.la

similarly you will need to have the lapack and the blas installed on
your system.

**********                                                **********
At this point you are ready to compile SIMPLE provided you have the
>>>> correct compilers installed on your system using if not got to
step (4) <<<<

$ perl Makefile_genAll.pl

It is recommneded that you read on to check that everything is setup
properly in particular for the GNU gcc and gfortran compilers suite.
**********                                                **********

(3) check that the accelerator variables are set to 0

$USE_MPI = 0;
$USE_GPU = 0;
$GPU_MODEL = 0;

and that all of the compilation directives are left as empty strings.

$DOPENMP = "-fopenmp";
$DCUDA = "";
$DOPENCL = "";
$DMAGMA = "";
$DMKL = "";
$DSIMPLE_MPI = "";
$DBENCH = "";

(4) if not already, specify the paths of the GNU compilers, make sure
that your are linking to the [GNU gcc] and **not the APPLE gcc**. On a
MacOSX systems it is typically found in:

$CC_COMPILER = "/usr/local/bin/gcc-4.9";
$GCC_COMPILER = "/usr/local/bin/g++-4.9";
$MPI_F_COMPILER = "/usr/loacl/bin/mpif90";
$FCOMPILER = "/usr/local/bin/gfortran-4.9";

if the GNU compiler suites are not installed already, install the
gfortran and gcc compiler from the GNU tool chain by following the
instruction from

https://wiki.helsinki.fi/display/HUGG/Installing+the+GNU+compilers+on+Mac+OS+X

and from

http://hpc.sourceforge.net/

select the bin tar ball 

gcc-5.0-bin.tar.gz, gfortran-5.0-bin.tar.gz (gfortran only), updated
Nov 2014 (Yosemite). 

and install it, if you do not have admin rights on your system then
point towards the installation folder.

After doing a 

$sudo make install

The gfortran and gcc compilers should appear in /usr/local/bin/

see that it works properly and check the version of the compiler

>$ /usr/local/bin/gfortran --version 
GNU Fortran (GCC) 5.0.0 20141005 (experimental)
Copyright (C) 2014 Free Software Foundation, Inc.

GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
You may redistribute copies of GNU Fortran
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING

similarly for the gcc compiler.

>$ /usr/local/bin/gcc --version
gcc (GCC) 5.0.0 20141005 (experimental)
Copyright (C) 2014 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

******
Be aware that Apple has also a gcc compiler which is different
than the GNU one and will generate errors at compilation.

>$ gcc --version 
Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/include/c++/4.2.1
Apple LLVM version 6.0 (clang-600.0.54) (based on LLVM 3.5svn)
Target: x86_64-apple-darwin13.4.0
Thread model: posix

you are pointing to the wrong compiler and you will get errors.

****** 

(5) Choose if you want to do a full compilation via

$ICOMPILE = {choose between 0,..,3, 3 is not yet finished};

0: will do make cleanall and a make and install everything from fresh
1: will only do a make and compile and link productions codes
2: will only generate the scripts and compile and link the production
codes 
3: same as 2: and is not yet implemented

(6) Once you have done all of that excute the PERL script via

$ perl Makefile_genAll.pl

------------
CPU and GPU
------------

(B) If you use CPU and GPU

first do a 

sudo chmod 777 /opt

so that you can create rwx sub-directories

create a folder make a 

/opt/Devel_tools

where you can install your 

NVIDIA_GPU_Computing_SDK
MAGMA-1.6.1 

tool kit and library repectively.

(1) do the installation as in (A) for the CPU part.

(2) include the CUDA and MAGMA code into the compilation by setting
the compilation directives in simple_user_input.pl script and set them
to

$DOPENMP = "-DOPENMP";         #Sets the OpenMP directives (default)
$DCUDA = "-DCUDA";             #Sets the CUDA driver
$DOPENCL = "-DOPENCL";         #Sets OpenCL not compatible with CUDA
$DMAGMA = "-DMAGMA";           #Sets the magma library and must go with CUDA
$DMKL = "-DMKL";               #Sets if we use MKL optimises lapacks
$DSIMPLE_MPI = "-DSIMPLE_MPI"; #this sets wether or not use MPI
$DBENCH = "-DBENCH";           #Sets the benchmarking flag

If you select -DCUDA no point selecting -DOPENCL

Set 

$USE_MPI = 0 or 1; #MPI not yet implemented
$USE_GPU = 1;

and choose your GPU model either the 0: CUDA or 1:OpenCL 

$GPU_MODEL = 0;

At the moment only CUDA is being imoplemented, OpenCL will come in
later release.

(2) If not already installed, install CUDA preferably CUDA-7.0 or
CUDA-7.5. Check that you have the correct hardware on your computer
and that it is compatible with your Operating System ie your MacOSX version

First, install the CUDA-toolkit from 

https://developer.nvidia.com/cuda-toolkit

which should install itself in directory 

/usr/local/cuda-7.x

We support the CUDA-7-0 and 7.5 version. And the path has a symbolic
link to  

/usr/local/cuda 

should be inserted in the

$CUDADIR="/usr/local/cuda";

You will need to install the NVIDIA_GPU_Computing_SDK from

http://developer.nvidia.com/cuda-toolkit-40 (ie http://developer.download.nvidia.com/compute/cuda/4_0/sdk/gpucomputingsdk_4.0.17_linux.run)

and for this you need to register as a developer. Pick "GPU Computing
SDK - complete package including all code samples". Check that the
driver has been installed properly from

modify your environment variables to 

#CUDA-7.0 Path to link to library
export CUDA_HOME=/usr/local/cuda-7.0
export PATH=$PATH:$CUDA_HOME/bin
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CUDA_HOME/lib

that is if you have upgraded from an earlier version otherwise 
the path is just /usr/local/cuda

You will then need to modify in simple_user_input.pm

$OPT_DEV_TOOLS_NVIDIA = "'/Developer/GPU Computing'";

Do not forget the extra ' ' otherwise it will not work as the space
between the U and C will not be detected and the path wil be broken.

cuDNN compatible with the cuda realease must be install into the
/usr/local/cuda7.x/lib64 for the *.so *.a and other and
/usr/local/cuda7.x/include for the cudnn.h file as a super user.

(3) if not already installed, install the Open MPI: Version 1.8.4

http://www.open-mpi.org/software/ompi/v1.8/

Copy the archive file to a temporary location to compile it. Open
Terminal.app and change directories to the new location of the Open
MPI archive. If you downloaded the source code as a .tar.gz file, you
can untar/zip it by typing 

cp {download folder}/openmpi-1.8.4.tar.gz /opt/Devel_tools

tar zxvf openmpi-1.8.4.tar.gz

which will create a directory called openmpi-1.8.4 

cd openmpi-1.8.4

You can now run the configuration script. If you only have one set of
compilers installed, you can run the configuration script by typing

./configure --prefix=/usr/local

If you have multiple compilers installed, you can specify which
compilers you would like to use as follows

./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=/usr/local

where you specify the C (CC), C++ (CXX), Fortran 77 (F77) and Fortran
90 and later (FC) compilers with the listed variables. Assuming the
configuration script runs without issue, you can compile the code by
typing make all which will compile the Open MPI libraries/binaries and
configure the wrappers for using the specified compilers. This should
take a bit... Again, assuming the code compiles without issue, you can
install Open MPI by typing 

sudo make install

Beware that using sudo can do major damage to your computer if you
aren't careful. To uninstall type make uninstall that will uninstall
the library.

(4) Install the MAGMA-1.6.1 library. Install the latest MAGMA library
from 

http://icl.cs.utk.edu/magma/software/index.html. 

Put it in 

/opt/Devel_tools/

then modify the variable in simple_user_input.pm file as

$MAGMADIR="/opt/Devel_tools/magma-1.6.1";

(5) open the simple_user_input.pm perl module file and edit the fields
that suits your installation there are the default usual values
already in the file but it allows alternative installation
configurations. A typical template file for MacOSX can be found in

./scripts/Template_MacOSX_10.9.5x64_CPU_simple_user_input.pm

the fields that need to be modified are and need to set properly for
the drivers and the compilers.

$SIMPLE_PATH= {path where simple is installed}
$ICOMPILE = {choose between 0,..,3} ;
# the name compiler default(linux): cc, gcc, gfortran, mpif90
# the name compiler default(MacOSX): /usr/local/bin/gcc, /usr/local/bin/gcc,
#                                    /usr/local/bin/gfortran, /usr/local/bin/mpif90
# The MPI directory default: MPI_DIR 
#[linux]: /usr        (Include): /usr/include/mpi
#[MacOSX]: /usr/local (include):/usr/local/include
$CC_COMPILER = {path to the cc compiler on MacOSX, default: /usr/local/bin/gcc}
$GCC_COMPILER = {path to the cc compiler on MacOSX, default: /usr/local/bin/g++}
$MPI_F_COMPILER = {path to the mpif90 compiler on MacOSX, default:/usr/local/bin/mpif90}
$FCOMPILER = {path to the gfortran or ifort compiler on MacOSX, default: /usr/local/bin/gfortran}

$MPI_DIR={path where the mpi directory is, default: /usr/local}
$MPI_DIR_INCLUDE={path where the mpi lib is installed, default: /usr/local/include"}
$USE_GPU = {wether or not GPU, default: 0}
$GPU_MODEL ={GPU model ONLY CUDA at the moment: default: 0}
$CUDADIR={path to the cuda driver, default: /usr/local/cuda}
$MAGMADIR={path to the MAGMA libL: default: /opt/Devel_tools/magma-1.6.1}}
$OPT_DEV_TOOLS_NVIDIA = {path to NVIDIA_GPU_SDK for OpenCL must have
the '' inside the "" otherwise it will not work!!!, default: '/Developer/GPU Computing'}
$OPENCL_SHARE_LIB ={path where lib on darwin is, if defaul;t install
then same as default, default: shared/lib/darwin";

$FFTW_LIB={for cluster need extra path after to the FFTW lib, default: "/usr/local/fftw/3.3.4-gcc/include/}

$OBJDIR={path to the object compiled files, default: obj/GFORTgpu}
$MODDIR={path to the .mod files, default: obj/GFORTgpu}

(6) Installing the cuDNN is very simple, login into your cuda account
and download the library from the

https://developer.nvidia.com/rdp/cudnn-download

if for cuda-7.0 cuDNN-v4 suffices, otherwise for cuda-7.5 cuDNN-v5 is
required. It may be instructive to get the samples as well.

Untar packages in a directory of your choice and a cuda directory will
appear. A file move is then required from both the lib64 and include
folders using

$sudo {mv or cp} include/cudnn.h /usr/local/cuda-7.0/include
$sudo {mv or cp} lib64/*.* /usr/local/cuda-7.0/lib64

similarly for the sample codes.

----------------------------------
Instalation (Linux)
-Ubuntu 14.10(cuda-7.0) or Ubuntu 15.04(cuda-7.5)
----------------------------------

(A) the installation is similar to the MacOSX above but the paths are
not the same.

Template files can be found in script for Ubuntu, Mint and SUSE

(1) Install CUDA. First, install the CUDA-toolkit from 

https://developer.nvidia.com/cuda-toolkit 

which should install itself in directory 

/usr/local/cuda 

Make sure to pick release 7 or newer. You will need to install the
NVIDIA_GPU_Computing_SDK from

http://developer.nvidia.com/cuda-toolkit-40

and for this you need to register as a developer. Pick "GPU Computing
SDK - complete package including all code samples".

#CUDA Path to link to library
export CUDA_HOME=/usr/local/cuda
export PATH=$PATH:$CUDA_HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib64

(2) Install OpenMPI. Install open mpi-2 or above 

(3) Install the Magma library. Install the latest MAGMA-1.6.1 library from 

http://icl.cs.utk.edu/magma/software/index.html. 

Put it in 

/opt/Devel_tools/

(4) Install the cuDNN library following the instruction above.

----------------------------------
Supported compilers
----------------------------------

(1) Compilers: 

GNU tool chain compilers

CPU: GNU 4.9 and above
GPU: GNU 4.9

are the supported compilers at the moment along as the NVIDIA nvcc compiler.

(2) anticipated compilers

Intel: ifort compiler
