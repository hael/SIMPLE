#!/usr/bin/perl
################################################################################
# Input variables for the compilation driving perl script                      #
#                                                                              #
#       Package: simple_user_input.pm                                          #
#                                                                              #
# This is the only place where the user needs to intervene in the compilation  #
# step of the library.                                                         #
################################################################################

package simple_user_input; # saved as simple_user_input.pm
use lib './scripts/';
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
use Config;
@ISA = ('Exporter');
@EXPORT = qw($ICOMPILE $SIMPLE_PATH $FCOMPILER $PLATFORM $CUDADIR $FFTW_LIB $FFTW_INC $OBJDIR $MODDIR $DEBUG $DEBUG_LEVEL $CC_COMPILER $GCC_COMPILER $DOPENMP $DCUDA $DBENCH $SET_OPTIMIZATION $op_sys $architecture getPlatform );

#####################################################################
# User-defined variables controlling compilation                    #
#####################################################################

# enter the SIMPLE root path
our$SIMPLE_PATH="/Users/creboul/Simple3";
# specifying the compiling directives
# with OpenMP:  -fopenmp, CUDA: -DCUDA
# Benchmarking: -DBENCH
our$DOPENMP = "-fopenmp";
our$DCUDA = "";
our$DBENCH = "";
# getting the platform details
our$PLATFORM = 1; # assuming Linux 
getPlatform();
# would you like to compile the library from scracth?
# make cleanall before make and link: 0   make only and generate scripts : 2
# make clean before make and link   : 1   link, compile production only  : 3 
our$ICOMPILE = 1;
# the name compiler default(linux): cc, gcc, gfortran
# the name compiler default(MacOSX): /usr/local/bin/gcc, /usr/local/bin/gcc,
#                                    /usr/local/bin/gfortran
our$CC_COMPILER = "";
our$GCC_COMPILER = "";
our$FCOMPILER = "gfortran";
# enter the CUDA libary path default: /usr/local/cuda
our$CUDADIR="";
# enter the fftw lib default: /usr/lib/x86_64-linux-gnu for [linux]
#                             /usr/local/lib for [MacOSX]
our$FFTW_LIB="/opt/local/lib/";
# on clusters we need extra path after module load fftw/3.3.4-gcc to identify the FFTW header
our$FFTW_INC="/opt/local/include/";
# Modules and objects directories. default: obj/SIMPLEOFILES
our$OBJDIR="obj/SIMPLEOFILES";
our$MODDIR="obj/SIMPLEOFILES";
# Set the optimization level: 0 1 2 or 3.
our$SET_OPTIMIZATION = 3;
# debug mode
# no debug mode: no, dubug mode: yes default: no
# if debug = yes the opmization level = null
our$DEBUG = "yes";
# debugging level "low" or "high"
our$DEBUG_LEVEL = "low";

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
