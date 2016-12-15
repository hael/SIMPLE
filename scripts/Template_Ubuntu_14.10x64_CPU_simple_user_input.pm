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
package simple_user_input; # saved as simple_user_input.pm
use lib './scripts/';
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
use Config;
use simple_clusterSpecs;

@ISA = ('Exporter');
@EXPORT = qw($ICOMPILE $SIMPLE_PATH $FCOMPILER $PLATFORM $USE_GPU $GPU_MODEL $CUDADIR $MAGMADIR $OPT_DEV_TOOLS_NVIDIA $OPENCL_SHARE_LIB $FFTW_LIB $FFTW_INC $OBJDIR $MODDIR $DEBUG $DEBUG_LEVEL $MPI_DIR $MPI_DIR_INCLUDE $CC_COMPILER $GCC_COMPILER $MPI_F_COMPILER $DOPENMP $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI $DBENCH $SET_OPTIMIZATION $op_sys $architecture getPlatform $SIMPLESYS %DISTR_ENV $EMAIL $NTHR $MEMSTR $SUBMITCMD $TIME_PER_IMAGE $SIMPLEBIN);

#####################################################################
# USER-DEFINED VARIABLES THAT CONTROL COMPILATION OF AND            #
# EXECUTION OF SIMPLE                                               #
#####################################################################
# enter the path of the application that you would like to run the SIMPLE from
our$SIMPLE_PATH="/mysimplepath/";
#specifying the compiling directives
#with OpenMP:  -fopenmp -DOPENMP,  CUDA: -DCUDA, MAGMA: -DMAGMA,
#     OpenCL: -DOPENCL,    Benchmarking: -DBENCH,
#        MKL: -DMKL,                MPI: -DSIMPLE_MPI
our$DOPENMP = "-fopenmp";
our$DCUDA = "";
our$DOPENCL = "";
our$DMAGMA = "";
our$DMKL = "";
our$DSIMPLE_MPI = "";
our$DBENCH = "";
#getting the platform details
our$PLATFORM = 1; # assuming Linux 
getPlatform();
# would you like to compile the library from scracth?
# make cleanall before make and link: 0   make only and generate scripts : 2
# make clean before make and link   : 1   link, compile production only  : 3 
our$ICOMPILE = 0;
# the name compiler default(linux): cc, gcc, gfortran, mpif90
# the name compiler default(MacOSX): /usr/local/bin/gcc, /usr/local/bin/gcc,
#                                    /usr/local/bin/gfortran, /usr/local/bin/mpif90
# The MPI directory default: MPI_DIR 
#[linux]: /usr        (include): /usr/include/mpi
#[MacOSX]: /usr/local (include):/usr/local/include
our$CC_COMPILER = "gcc";
our$GCC_COMPILER = "g++";
our$MPI_F_COMPILER = "/usr/bin/mpif90";
our$FCOMPILER = "gfortran";
#our$MPI_DIR="/usr";
#our$MPI_DIR_INCLUDE="/usr/include/mpi";
our$MPI_DIR="";
our$MPI_DIR_INCLUDE="";
#enter wether or not you want to use MPI compilation
# 0: No or 1: Yes
our$USE_MPI = 0;
#enter the compilation model that you would like to use
# 0: 100% CPU or 1: Hybrid CPU/GPU
our$USE_GPU = 0;
# Which GPU modelling sets the driver tolink the GPU code
# either the 0: CUDA or 1:OpenCL 
our$GPU_MODEL = 0;
#enter the CUDA libary path default: /usr/local/cuda
#our$CUDADIR="/usr/local/cuda";
our$CUDADIR="";
#The Magma directory default:/opt/Devel_tools/magma-1.6.1
#our$MAGMADIR="/opt/Devel_tools/magma-1.6.1";
our$MAGMADIR="";
# The OpenCL paths for the SDK#####
#for [linux]
#our$OPT_DEV_TOOLS_NVIDIA = "/opt/Devel_tools/NVIDIA_GPU_Computing_SDK";
our$OPT_DEV_TOOLS_NVIDIA = "";
#for [MacOSX]
#$OPT_DEV_TOOLS_NVIDIA = "'/Developer/GPU Computing'";
# The lib to be used is set here default: shared/lib/linux 
# $OPT_DEV_TOOLS_NVIDIA and are concatenated together for linking the library 
#for [linux]
#our$OPENCL_SHARE_LIB = "shared/lib/linux";
our$OPENCL_SHARE_LIB = "";
#for [MacOSX]
#$OPENCL_SHARE_LIB = "shared/lib/darwin";
#The fftw lib default: /usr/lib/x86_64-linux-gnu
#for [linux]
our$FFTW_LIB="/usr/lib/x86_64-linux-gnu";
#for [MacOSX]
#$FFTW_LIB="/usr/local/lib";
#for cluster need extra path after 
# > module load fftw/3.3.4-gcc ;
our$FFTW_INC="/usr/local/fftw/3.3.4-gcc/include/"; #[Massive]
# Modules and objects directories. default: obj/GFORTgpu
our$OBJDIR="obj/GFORTgpu";
our$MODDIR="obj/GFORTgpu";
#Set the optimization level: 0 1 2 or 3.
our$SET_OPTIMIZATION = 3;
#debug mode
# no debug mode: "no", dubug mode: "yes" default: no
# if debug = yes the opmization level = null
our$DEBUG = "no";
#debugging level "low" or "high"
our$DEBUG_LEVEL = "low";

#####################################################################
# USER-DEFINED VARIABLES THAT CONTROL WHICH CLUSTER ENVIRONMENT TO  #
# USE AND HOW TO USE THE RESOURCES                                  #
#####################################################################
our$SIMPLESYS      = 'LOCAL';            # Name of system
our%DISTR_ENV       = %LOCAL_DISTR_ENV;  # Defines the environment for distributed execution     
our$EMAIL           = 'myemail@uni.edu'; # e-mail for failure report
our$NTHR            = 4;                 # number of threads (CPUs per core)
our$MEMSTR          = '500';             # string descriptor for memory
our$TIME_PER_IMAGE = 100;                # time per image (in seconds)

# SETTINGS FOR COMLIN_CLUSTER ON MASSIVE2
# $NTHR           = 1
# $MEMSTR         = '500';
# $TIME_PER_IMAGE = 5;

# SETTINGS FOR PRIME2/PRIME2_CLUSTER ON MASSIVE2
# $NTHR           = 8
# $MEMSTR         = '32000';
# $TIME_PER_IMAGE = 100;

# SETTINGS FOR PRIME2/PRIME2_CLUSTER IN OXFORD
# $NTHR           = 8
# $MEMSTR         = '9gb';
# $TIME_PER_IMAGE = =100;

#####################################################################
# OTHER VARIABLES TO EXPORT                                         #
#####################################################################
our$SIMPLEBIN = $SIMPLE_PATH.'/bin'; # location of compiled binaries
our$SUBMITCMD = $DISTR_ENV{'SUBMITCMD'};   

#####################################################################
# Subroutine to get the platform details.                           #
#####################################################################
sub getPlatform {
    #extracting the Platrform details
    our$op_sys = $Config{osname};
    our$architecture = $Config{archname};
    #The platform gets detected from the system using that you wish to run on
    # 0: MacOsX    1: Linux
    if( $op_sys eq 'darwin' ){
	    $PLATFORM = 0;
    }elsif( $op_sys eq 'linux' ){
	    $PLATFORM = 1;
    }
}
1;
