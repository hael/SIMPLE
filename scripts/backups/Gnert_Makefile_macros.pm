#!/usr/bin/perl
################################################################################
# Script that generates the Makefile_macros for a given achitecture            #
# and Input variables for the compilation driving perl script from the         #
#                                                                              #
# use package simple_user_input                                                #
#                                                                              #
#       Package Gnert_Makefile_macros                                          #
#                                                                              #
#perl script that lets the user input the global variables that defines the    #
#global variables for the entire SIMPLE Application. This is the only place    #
#where the user needs to intervene in the compilation step of the library.     #
################################################################################

package Gnert_Makefile_macros;  # saved as simple_user_input.pm
use lib './';
use Exporter;
use warnings;
use Config;
################################################################################
# Calling the input module file simple_user_input.pm                           #
#   **** Do not modify this script, this sript automates the compilation ****  #
#   ****                         procedure                               ****  #
#                                                                              #
# The input variables are to be modifed in the                                 #
#                                                                              #
#                          simple_user_input.pm                                #
#                                                                              #
################################################################################

use simple_user_input;

################################################################################

@ISA = ('Exporter');
@EXPORT = qw(make_Makefile_macros);

#Opening the files for the 
$filename_mkma = 'temp_gen_Makefile_macros';
open($mkma, '>', $filename) or die "Could not open file '$filename' $!";

make_Makefile_macros();

close $mkma;
################################################################################
# Soubroutines for file handling.                                              #
#                                                                              #
################################################################################
#
sub make_Makefile_macros {
    #($option_Makefile_macros_f90_flags, $option_Makefile_macros_f77_flags, $option_Makefile_macros_mpif90_flags, $option_Makefile_macros_mpif77_flags) = @_;

    print $mkma qq[# Makefile_macros for the g95 fortran compiler.\n];
    print $mkma qq[\n];
    print $mkma qq[####### The project path \#####\n];
    print $mkma qq[\n];
    print $mkma qq[Simple_source=$SIMPLE_PATH\n];
    print $mkma qq[\n];
    print $mkma qq[####### The switches \#####\n];
    print $mkma qq[\n];
    print $mkma qq[DLIB = -DCUDA\n];
    print $mkma qq[#DLIB = -DMKL\n];
    print $mkma qq[DMAGMA = -DMAGMA\n];
    print $mkma qq[DSIMPLE_MPI = -DSIMPLE_MPI\n];
    print $mkma qq[\n];
    print $mkma qq[####### The MPI paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[MPIDIR=$MPI_DIR\n];
    print $mkma qq[MPIDIR_INCLUDE=\$(MPIDIR)/include/mpi\n];
    print $mkma qq[\n];
    print $mkma qq[####### The compilers paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CCCOMP=$CC_COMPILER\n];
    print $mkma qq[GCC=$GCC_COMPILER\n];
    print $mkma qq[MPIFORTRAN=$MPI_F_COMPILER\n];
    print $mkma qq[GFORTRAN=$FCOMPILER\n];
    print $mkma qq[\n];
    print $mkma qq[####### The CUDA and MAGMA paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CUDADIR=$CUDADIR\n];
    print $mkma qq[\n];
    print $mkma qq[MAGMADIR=$MAGMADIR\n];
    print $mkma qq[MAGMADIR_CONTROL=\$(MAGMADIR)/control\n];
    print $mkma qq[MAGMADIR_INCLUDE=\$(MAGMADIR)/include\n];
    print $mkma qq[\n];
    print $mkma qq[####### The OpenCL paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[OPT_DEV_TOOLS_NVIDIA = $OPT_DEV_TOOLS_NVIDIA\n];
    print $mkma qq[\n];
    print $mkma qq[####### Set MAGMA-ADDONS={Magma version} to compile Magma addons. \#####\n];
    print $mkma qq[\n];
    print $mkma qq[GPU_MODEL=0.0     \#0: tesla arch, 1: fermi arch \n];
    print $mkma qq[\n];
    print $mkma qq[####### Modules and objects directories. \#####\n];
    print $mkma qq[\n];
    print $mkma qq[OBJDIR=$OBJDIR\n];
    print $mkma qq[MODDIR=$MODDIR\n];
    print $mkma qq[\n];
    print $mkma qq[##############\n];
    print $mkma qq[# C compiler.#\n];
    print $mkma qq[##############\n];
    print $mkma qq[\n];
    print $mkma qq[CC=\$(CCCOMP)\n];
    print $mkma qq[CFLAGS= -O3 -O -DADD_ -DBENCH \$(DLIB)\n];
    print $mkma qq[CCLIB=\$(CC) -c \$(CFLAGS) -I \$(MODDIR)                                            \\\n];
    print $mkma qq[                         -I \$(OPT_DEV_TOOLS_NVIDIA)/OpenCL/common/inc/           \\\n];
    print $mkma qq[                         -I \$(OPT_DEV_TOOLS_NVIDIA)/shared/inc/                  \\\n];
    print $mkma qq[                         -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                         -I \$(CUDADIR)/src                                       \\\n];
    print $mkma qq[                         -I \$(MAGMADIR_INCLUDE)                                  \\\n];
    print $mkma qq[                         -I \$(MAGMADIR_CONTROL)                                  \\\n];
    print $mkma qq[                         -I \$(MPIDIR_INCLUDE)                                    \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include                             \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/simple                      \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/OpenCL                      \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/mpe\n];
    print $mkma qq[\n];
    print $mkma qq[################\n];
    print $mkma qq[# C++ compiler.#\n];
    print $mkma qq[################\n];
    print $mkma qq[\n];
    print $mkma qq[CPP=\$(GCC)\n];
    print $mkma qq[CPPFLAGS= -O3 -DADD_ \$(DMAGMA) \$(DLIB)\n];
    print $mkma qq[CPPCLIB=\$(CPP) -c \$(CPPFLAGS) -I \$(MODDIR)                                            \\\n];
    print $mkma qq[                              -I \$(OPT_DEV_TOOLS_NVIDIA)/OpenCL/common/inc/           \\\n];
    print $mkma qq[                              -I \$(OPT_DEV_TOOLS_NVIDIA)/shared/inc/                  \\\n];
    print $mkma qq[                              -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                              -I \$(CUDADIR)/src                                       \\\n];
    print $mkma qq[                              -I \$(MAGMADIR_INCLUDE)                                  \\\n];
    print $mkma qq[                              -I \$(MAGMADIR_CONTROL)                                  \\\n];
    print $mkma qq[                              -I \$(MPIDIR_INCLUDE)                                    \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include                             \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/simple                      \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/OpenCL                      \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/mpe                         \\\n];
    print $mkma qq[                              -I \$(Simple_source)/testCode\n];
    print $mkma qq[\n];
    print $mkma qq[#################\n];
    print $mkma qq[# CUDA compiler.#\n];
    print $mkma qq[#################\n];
    print $mkma qq[\n];
    print $mkma qq[#ifeq (\$(GPU_MODEL),0.0)\n];
    print $mkma qq[#	PUSHMEM_GPU= -arch sm_13 -DGPUSHMEM=130 -gencode arch=compute_13,code=compute_13 -gencode arch=compute_10,code=compute_10\n];
    print $mkma qq[#else\n];
    print $mkma qq[        PUSHMEM_GPU= -arch sm_20 -DGPUSHMEM=200 -gencode arch=compute_20,code=compute_20\n];
    print $mkma qq[#endif\n];
    print $mkma qq[\n];
    print $mkma qq[NVCC=\$(CUDADIR)/bin/nvcc\n];
    print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3\n];
    print $mkma qq[NVCCCLIB=\$(NVCC) -c \$(NVCCFLAGS) -I \$(MODDIR)                               \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include/cuda           \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include                \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include/simple         \\\n];
    print $mkma qq[                                 -I \$(CUDADIR)/include                      \\\n];
    print $mkma qq[                                 -I \$(MAGMADIR_INCLUDE)                     \\\n];
    print $mkma qq[                                 -I \$(MAGMADIR_CONTROL)\n];
    print $mkma qq[\n];
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# Fortran 77 compiler.#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[F77C=\$(GFORTRAN)\n];
    print $mkma qq[F77FLAGS=-O3 -cpp\n];
    print $mkma qq[F77CLIB=\$(F77C) -c \$(F77FLAGS)\n];
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# Fortran 90 compiler.#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[# \$(DLIB)\n];
    print $mkma qq[F90C=\$(GFORTRAN)\n];
    print $mkma qq[F90FLAGS= $option_Makefile_macros_f90_flags\n];
    print $mkma qq[F90FLAGS77=$option_Makefile_macros_f77_flags\n];
    print $mkma qq[F90CLIB=\$(F90C) -c \$(F90FLAGS) -I .                                \\\n];
    print $mkma qq[                               -I \$(MODDIR)                        \\\n];
    print $mkma qq[                               -J \$(MODDIR)                        \\\n];
    print $mkma qq[                               -I \$(Simple_source)/include/simple  \\\n];
    print $mkma qq[                               -I \$(MPIDIR_INCLUDE)                \\\n];
    print $mkma qq[                               -I \$(MAGMADIR_INCLUDE)              \\\n];
    print $mkma qq[                               -I \$(MAGMADIR_CONTROL)\n];
    print $mkma qq[\n];
    print $mkma qq[\n];
    print $mkma qq[F90CLIB77=\$(F90C) -c \$(F90FLAGS77) -I .                                \\\n];
    print $mkma qq[                                   -I \$(MODDIR)                        \\\n];
    print $mkma qq[                                   -J \$(MODDIR)                        \\\n];
    print $mkma qq[                                   -I \$(Simple_source)/include/simple  \\\n];
    print $mkma qq[                                   -I \$(MPIDIR_INCLUDE)                \\\n];
    print $mkma qq[                                   -I \$(MAGMADIR_INCLUDE)              \\\n];
    print $mkma qq[                                   -I \$(MAGMADIR_CONTROL)\n];
    print $mkma qq[\n];
    print $mkma qq[F90POST=\n];
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# MPIF90 compiler.    \#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90C=\$(MPIFORTRAN)\n];
    print $mkma qq[MPIF90FLAGS=$option_Makefile_macros_mpif90_flags\n];
    print $mkma qq[MPIF90FLAGS77=$option_Makefile_macros_mpif77_flags\n];
    print $mkma qq[MPIF90CLIB=\$(MPIF90C) -c \$(MPIF90FLAGS) -I .                      \\\n];
    print $mkma qq[                               -I \$(MODDIR)                       \\\n];
    print $mkma qq[                               -J \$(MODDIR)                       \\\n];
    print $mkma qq[                               -I \$(MPIDIR_INCLUDE)               \\\n];
    print $mkma qq[                               -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90CLIB77=\$(MPIF90C) -c \$(MPIF90FLAGS77) -I .                     \\\n];
    print $mkma qq[                                  -I \$(MODDIR)                       \\\n];
    print $mkma qq[                                  -J \$(MODDIR)                       \\\n];
    print $mkma qq[                                  -I \$(MPIDIR_INCLUDE)               \\\n];
    print $mkma qq[                                  -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90POST=\n];
    print $mkma qq[\n];
}
1;

