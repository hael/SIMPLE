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
@EXPORT = qw($ICOMPILE $SIMPLE_PATH $FCOMPILER $PLATFORM $USE_GPU $GPU_MODEL $CUDADIR $MAGMADIR $OPT_DEV_TOOLS_NVIDIA $OPENCL_SHARE_LIB $FFTW_LIB $FFTW_INC $OBJDIR $MODDIR $DEBUG $DEBUG_LEVEL $MPI_DIR $MPI_DIR_INCLUDE $CC_COMPILER $GCC_COMPILER $MPI_F_COMPILER $DOPENMP $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI $DBENCH $SET_OPTIMIZATION $op_sys $architecture getPlatform);

#enter the path of the application that you would like to run the SIMPLE from
$SIMPLE_PATH="/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet";

#specifying the compiling directives
#with OpenMP:  -fopenmp -DOPENMP,  CUDA: -DCUDA, MAGMA: -DMAGMA,
#     OpenCL: -DOPENCL,    Benchmarking: -DBENCH,
#        MKL: -DMKL,                MPI: -DSIMPLE_MPI
$DOPENMP = "-fopenmp -DOPENMP";
$DCUDA = "-DCUDA";
$DOPENCL = "-DOPENCL";
$DMAGMA = "-DMAGMA";
$DMKL = "-DMKL";
$DSIMPLE_MPI = "-DSIMPLE_MPI";
$DBENCH = "-DBENCH";

#getting the platform details
getPlatform();

# would you like to compile the library from scracth?
# make cleanall before make and link: 0   make only and generate scripts : 2
# make clean before make and link   : 1   link, compile production only  : 3 
$ICOMPILE = 0;

# the name compiler default(linux): cc, gcc, gfortran, mpif90
# the name compiler default(MacOSX): /usr/local/bin/gcc, /usr/local/bin/gcc,
#                                    /usr/local/bin/gfortran, /usr/local/bin/mpif90
# The MPI directory default: MPI_DIR 
#[linux]: /usr        (include): /usr/include/mpi
#[MacOSX]: /usr/local (include):/usr/local/include

$CC_COMPILER = "gcc";
$GCC_COMPILER = "g++";
$MPI_F_COMPILER = "/usr/bin/mpif90";

$FCOMPILER = "gfortran";

$MPI_DIR="/usr";
$MPI_DIR_INCLUDE="/usr/include/mpi";
#$MPI_DIR="";
#$MPI_DIR_INCLUDE="";

#enter wether or not you want to use MPI compilation
# 0: No or 1: Yes
$USE_MPI = 0;
#enter the compilation model that you would like to use
# 0: 100% CPU or 1: Hybrid CPU/GPU
$USE_GPU = 1;
# Which GPU modelling sets the driver tolink the GPU code
# either the 0: CUDA or 1:OpenCL 
$GPU_MODEL = 0;

#enter the CUDA libary path default: /usr/local/cuda
$CUDADIR="/usr/local/cuda";
#$CUDADIR="";

#The Magma directory default:/opt/Devel_tools/magma-1.6.1
$MAGMADIR="/opt/Devel_tools/magma-1.6.1";
#$MAGMADIR="";

# The OpenCL paths for the SDK#####
#for [linux]
$OPT_DEV_TOOLS_NVIDIA = "/opt/Devel_tools/NVIDIA_GPU_Computing_SDK";
#$OPT_DEV_TOOLS_NVIDIA = "";
#for [MacOSX]
#$OPT_DEV_TOOLS_NVIDIA = "'/Developer/GPU Computing'";

# The lib to be used is set here default: shared/lib/linux 
# $OPT_DEV_TOOLS_NVIDIA and are concatenated together for linking the library 
#for [linux]
$OPENCL_SHARE_LIB = "shared/lib/linux";
#$OPENCL_SHARE_LIB = "";
#for [MacOSX]
#$OPENCL_SHARE_LIB = "shared/lib/darwin";

#The fftw lib default: /usr/lib/x86_64-linux-gnu
#for [linux]
$FFTW_LIB="/usr/lib/x86_64-linux-gnu";
#for [MacOSX]
#$FFTW_LIB="/usr/local/lib";
#for cluster need extra path after 
# > module load fftw/3.3.4-gcc ;
$FFTW_INC="/usr/local/fftw/3.3.4-gcc/include/"; #[Massive]

# Modules and objects directories. default: obj/GFORTgpu
$OBJDIR="obj/GFORTgpu";
$MODDIR="obj/GFORTgpu";

#Set the optimization level: 0 1 2 or 3.
$SET_OPTIMIZATION = 3;

#debug mode
# no debug mode: "no", dubug mode: "yes" default: no
# if debug = yes the opmization level = null
$DEBUG = "no";
#debugging level "low" or "high"
$DEBUG_LEVEL = "low";

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

