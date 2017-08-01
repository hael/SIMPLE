#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Function to print header for the SIMPLE library for the bench.py
# Author: Frederic Bonnet
# Date: 23/07/2016
#----------------------------------------------------------------------------
#System tools
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
################################################################################
# The helper method use at the command line                                    #
################################################################################
def help():
    print "\033[0;94m **************** start check_help ****************************"
    print sys.argv
    print""
    print"\033[0;94m Usage: make [OPTIONS]... [VAR=VALUE]..."
    print""
    print"\033[0;94m [OPTIONS]:"
    print"\033[0;94m check            Runs ALL the checks 1WCM, Clustering, check_cpu"
    print"\033[0;94m bench            Runs a set of bench marking codes, can take a long time!!!"
    print"\033[0;94m check_help       Prints help options for the different options"
    print"\033[0;94m check_news       Prints the news about the release and upcoming release"
    print"\033[0;94m check_cpu        Runs ALL the checks cpu related, FFTW-1,2,3D, lapacks-Inversion"
    print"\033[0;94m check_gpu        Runs ALL the checks gpu related, cuFFT-1,2,3D, cublas, MAGMA Inversion"
    print"\033[0;94m wc               Word counts of the code"
    print"\033[0;94m tar              Creates a tar ball of the library"
    print""
    print"\033[0;94m clean            Cleans all the objet files in /obj/GFORTgpu"
    print"\033[0;94m cleanall         Cleans all the objet files in /obj/GFORTgpu and all compiled code"
    print"\033[0;94m clean_check_cpu  Cleans objet files from the make check_cpu"
    print"\033[0;94m clean_check_gpu  Cleans objet files from the make check_cpu"
    print""
    print"\033[0;94m [VAR]:"
    print"\033[0;94m use_gpu          Use GPU (CUDA), yes or no. Default no for PRIME3D"
    print"\033[0;94m bench_gpu        Use CPU to compute the equivalent to (CUDA), yes or no. Default no for PRIME3D"
    print"\033[0;94m fix_gpu          Specify if fixing onto a gpu, must be used with set_gpu, yes or no. Default no for PRIME3D"
    print"\033[0;94m set_gpu          Specify the device number used for comoputation, 0,..,MAX_N_GPU=8. Default 0 for PRIME3D"
    print"\033[0;94m help             Prints help options for the different options"
    print""
    print"\033[0;94m VALUE:"
    print"\033[0;94m yes/no"
    print"\033[0;94m {0,..,MAX_N_GPU}"
    print""
    return
