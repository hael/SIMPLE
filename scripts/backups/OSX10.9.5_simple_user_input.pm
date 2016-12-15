#!/usr/bin/perl
################################################################################
# Input variables for the compilation driving perl script                      #
#                                                                              #
#       Package simple_user_input                                              #
#                                                                              #
# perl script that lets the user input the global variables that defines the   #
# global variables for the entire SIMPLE Application. This is the only place   #
# where the user needs to intervene in the compilation step of the library.    #
################################################################################
#
package simple_user_input;  # saved as simple_user_input.pm
use lib './';
use Exporter;
use warnings;
use Config;

@ISA = ('Exporter');
@EXPORT = qw($ICOMPILE $SIMPLE_PATH $FCOMPILER $PLATFORM $USE_GPU $GPU_MODEL $CUDADIR $MAGMADIR $OPT_DEV_TOOLS_NVIDIA $OPENCL_SHARE_LIB $FFTW_LIB $OBJDIR $MODDIR $DEBUG $MPI_DIR $MPI_DIR_INCLUDE $CC_COMPILER $GCC_COMPILER $MPI_F_COMPILER $option_mkma_f90_flags $option_mkma_f77_flags $option_mkma_mpif90_flags $option_mkma_mpif77_flags $option_in $op_sys $architecture getPlatform);

#enter the path of the application that you would like to run the SIMPLE from
$SIMPLE_PATH="/Users/frederic/SourceCode/temp/Simple_Restruct.projet";
#would you like to compile the library from scracth?
# make cleanall before make: 0, make: 1, generate scripts only: 2, compile production only: 3 
$ICOMPILE = 0;

#getting the platform details
getPlatform();

# the name compiler default(linux): cc, gcc, gfortran, mpif90
# the name compiler default(MacOSX): /usr/local/bin/gcc, /usr/local/bin/gcc,
#                                    /usr/local/bin/gfortran, /usr/local/bin/mpif90
# The MPI directory default: MPI_DIR 
#[linux]: /usr        (include): /usr/include/mpi
#[MacOSX]: /usr/local (include):/usr/local/include

$CC_COMPILER = "/usr/local/bin/gcc";
$GCC_COMPILER = "/usr/local/bin/gcc";
$MPI_F_COMPILER = "/usr/local/bin/mpif90";

$FCOMPILER = "/usr/local/bin/gfortran";

$MPI_DIR="/usr/local";
$MPI_DIR_INCLUDE="/usr/local/include";

#enter the compilation model that you would like to use
# 0: 100% CPU or 1: Hybrid CPU/GPU
$USE_GPU = 0;
# Which GPU modelling sets the driver tolink the GPU code
# either the 0: CUDA or 1:OpenCL 
$GPU_MODEL = 0;

#enter the CUDA libary path default: /usr/local/cuda
$CUDADIR="/usr/local/cuda";

#The Magma directory default:/opt/Devel_tools/magma-1.6.1
$MAGMADIR="/opt/Devel_tools/magma-1.6.1";

# The OpenCL paths for the SDK#####
#for [linux]
#$OPT_DEV_TOOLS_NVIDIA = "/opt/Devel_tools/NVIDIA_GPU_Computing_SDK";
#for [MacOSX]
$OPT_DEV_TOOLS_NVIDIA = "'/Developer/GPU Computing'";
# The lib to be used is set here default: shared/lib/linux 
# $OPT_DEV_TOOLS_NVIDIA and are concatenated together for linking the library 
#for [linux]
#$OPENCL_SHARE_LIB = "shared/lib/linux";
#for [MacOSX]
$OPENCL_SHARE_LIB = "shared/lib/darwin";

#The fftw lib default: /usr/lib/x86_64-linux-gnu
#for [linux]
#$FFTW_LIB="/usr/lib/x86_64-linux-gnu";
#for [MacOSX]
$FFTW_LIB="/usr/local/lib";

# Modules and objects directories. default: obj/GFORTgpu
$OBJDIR="obj/GFORTgpu";
$MODDIR="obj/GFORTgpu";

#compiler optimisation varaibales (This goes into the Makefile_macros file) this options are for linux
#for [linux]
#linking options (this goes into the f90g95_local file)
#$option_in                = '-ffree-form -cpp -O3 -fPIC -fno-second-underscore -fopenmp -DBENCH -DMAGMA';
#compiling options for the Makefile_macros
#$option_mkma_f90_flags    = '-ffree-form -cpp -O3 -fno-second-underscore -DBENCH $(DMAGMA) $(DSIMPLE_MPI) $(DLIB)';
#$option_mkma_f77_flags    = '-ffixed-form -cpp -O3 -fPIC -fno-second-underscore -DBENCH $(DMAGMA) $(DSIMPLE_MPI) $(DLIB)';
#$option_mkma_mpif90_flags = '-ffree-form -g -cpp -O3 -fno-second-underscore -DBENCH';
#$option_mkma_mpif77_flags = '-ffixed-form -g -cpp -O3 -fPIC -fno-second-underscore';
#for [MacOSX]
$option_in                = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -O3 -fopenmp -DBENCH -DMAGMA';
#compiling options for the Makefile_macros
$option_mkma_f90_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -DBENCH $(DMAGMA) $(DSIMPLE_MPI) $(DLIB)';
$option_mkma_f77_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -DBENCH $(DMAGMA) $(DSIMPLE_MPI) $(DLIB)';
$option_mkma_mpif77_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
$option_mkma_mpif90_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fno-second-underscore';

#debug mode
# no debug mode: no, dubug mode: yes default: no
$DEBUG = "no";

################################################################################
# Subroutine to get the platform details.                                      #
#                                                                              #
################################################################################
# 
#putting the variable for the module export
sub getPlatform {

    #extracting the Platrform details
    $op_sys = $Config{osname};
    $architecture = $Config{archname};
    #The platform gets detected from the system using that you wish to run on
    # 0: MacOsX    1: Linux
    if( $op_sys eq 'darwin' ){
	$PLATFORM = 0;
    }elsif( $op_sys eq 'linux' ){
	$PLATFORM = 1;
    }

}

1;

