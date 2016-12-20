#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Build the final report (read *.report)
# Author: Frederic Bonnet
# Date: 18/01/2015
#----------------------------------------------------------------------------

import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
#from psutil import virtual_memory
#import yaml
def print_header():
    print "\033[0;94m *********************************************************"
    print "\033[0;94m *This is the python sript to run the SIMPLE check suite.*"
    print "\033[0;94m *                                                       *"
    print "\033[0;94m *  ####      #    #    #  #####   #       ######        *"
    print "\033[0;94m * #          #    ##  ##  #    #  #       #             *"
    print "\033[0;94m *  ####      #    # ## #  #    #  #       #####         *"
    print "\033[0;94m *      #     #    #    #  #####   #       #             *"
    print "\033[0;94m * #    #     #    #    #  #       #       #             *"
    print "\033[0;94m *  ####      #    #    #  #       ######  ######        *"
    print "\033[0;94m *                                                       *"
    print "\033[0;94m *********************************************************"
    print
    return
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
def myplatform():
    print "\033[0;94m **************** Platform details ****************************"
    print platform.version()
    print platform.platform()
    print platform.uname()
    #print platform.system()
    #print platform.processor()
    print ""
    cpu_count = multiprocessing.cpu_count()
    print "number of avaible processors: ",cpu_count
    return 8 #cpu_count

print_header()
cpu_count = myplatform()
len_sys = len(sys.argv)
if (len_sys == 2):
    if sys.argv[1] in ['--help']: help()
    exit(0)
if (len_sys > 2):
    if sys.argv[1] in ['--use_gpu=yes' , '--use_gpu=no']: print "sys.argv[1]: ",sys.argv[1]
    if sys.argv[2] in ['--use_gpu=yes' , '--use_gpu=no']: print "sys.argv[2]: ",sys.argv[2]
    if sys.argv[3] in ['--use_gpu=yes' , '--use_gpu=no']: print "sys.argv[3]: ",sys.argv[3]

    if sys.argv[1] in ['--bench_gpu=yes' , '--bench_gpu=no']: print "sys.argv[1]: ",sys.argv[1]
    if sys.argv[2] in ['--bench_gpu=yes' , '--bench_gpu=no']: print "sys.argv[2]: ",sys.argv[2]
    if sys.argv[3] in ['--bench_gpu=yes' , '--bench_gpu=no']: print "sys.argv[3]: ",sys.argv[3]

    if sys.argv[1] in ['--help=yes' , '--help=no']: print "sys.argv[1]: ",sys.argv[1]
    if sys.argv[2] in ['--help=yes' , '--help=no']: print "sys.argv[2]: ",sys.argv[2]
    if sys.argv[3] in ['--help=yes' , '--help=no']: print "sys.argv[3]: ",sys.argv[3]
    if sys.argv[1] in ['--help=yes']: help()
    if sys.argv[2] in ['--help=yes']: help()
    if sys.argv[3] in ['--help=yes']: help()

LOCAL_DIR = os.getcwd()
print "\033[0;92m ********************Starting the checks******************"

#if (len_sys == 3):
#!!!!!!!!!!!!
# 1WCM testcode
#!!!!!!!!!!!!
if ( sys.argv[1] not in ['--use_gpu=yes' , '--use_gpu=no', '--bench_gpu=yes' , '--bench_gpu=no', '--help=yes' , '--help=no'] and sys.argv[2] not in ['--use_gpu=yes' , '--use_gpu=no', '--bench_gpu=yes' , '--bench_gpu=no', '--help=yes' , '--help=no']):
        print sys.argv[1]
        print sys.argv[2]
        print sys.argv[3]
        print "\033[0;92m ******************** 1WCM... ****************************\033[0m"
        print "\033[0;92m ***** no CPU or GPU specifier defaulting on CPU *********\033[0m"
        os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 no no; rm polii.spi'%cpu_count)
else:
    print "\033[0;92m ******************** 1WCM... ****************************\033[0m"
    if ( sys.argv[1] in ['--use_gpu=no' , '--use_gpu='] ):
        if ( sys.argv[2] not in ['--bench_gpu=yes' , '--bench_gpu=no'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 no no; rm polii.spi'%cpu_count)
        if ( sys.argv[2] in ['--bench_gpu=yes'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 no yes; rm polii.spi'%cpu_count)
        if ( sys.argv[2] in ['--bench_gpu=no'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 no no; rm polii.spi'%cpu_count)
    if ( sys.argv[1] in ['--use_gpu=yes'] ):
        if ( sys.argv[2] not in ['--bench_gpu=yes' , '--bench_gpu=no'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 yes no; rm polii.spi'%cpu_count)
        if ( sys.argv[2] in ['--bench_gpu=yes'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            print "\033[0;91m *********************************************************\033[0m"
            print "\033[0;91m Exit: not coherent to run bench_gpu=yes with use_gpu=yes \033[0m"
            print "\033[0;91m *********************************************************\033[0m"
        if ( sys.argv[2] in ['--bench_gpu=no'] ):
            print sys.argv[1]
            print sys.argv[2]
            print sys.argv[3]
            print "sys.argv[1]: ",sys.argv[1],"sys.argv[2]: ",sys.argv[2]
            os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 32 5 %i 240 1.77 76 yes no; rm polii.spi'%cpu_count)
    
#!!!!!!!!!!!!
# GPU Clustering
#!!!!!!!!!!!!
print "\033[0;92m ******************* Clustering... ***********************\033[0m"
#os.system('cd ./checks/simple_LibTests/Clustering/ ; rm *.txt *.spi ; cp ../1WCM/spi.polii . ; ln -s spi.polii polii.spi ; csh create_simulData_clutering.csh ./polii.spi ./ 32 5 %i 240 1.77 76 ; rm polii.spi spi.polii'%cpu_count)
#!!!!!!!!!!!!
# GPU Symmetry
#!!!!!!!!!!!!
#    print "\033[0;92m ******************* Symmetry... *************************\033[0m"
#!!!!!!!!!!!!
# CPU testcode
#!!!!!!!!!!!!
print "\033[0;92m **************** Test_codes (CPU)... ********************\033[0m"
#os.system('cd ./checks/simple_LibTests/cpu_test_code ; make ; csh ./run_Testing_cpu.csh')
#!!!!!!!!!!!!
# GPU testcode
#!!!!!!!!!!!!
if sys.argv[1] in ['--use_gpu=yes']:
    print "\033[0;92m **************** Test_codes (GPU)... ********************\033[0m"
#    os.system('cd ./checks/simple_LibTests/gpu_test_code ; make ; csh ./run_Testing_gpu.csh')
elif sys.argv[2] in ['--use_gpu=yes']:
    print "\033[0;92m **************** Test_codes (GPU)... ********************\033[0m"
#    os.system('cd ./checks/simple_LibTests/gpu_test_code ; make ; csh ./run_Testing_gpu.csh')
    
#Check if the output is a tty to print in colour
if sys.stdout.isatty():
    start_fail = "\033[0;31m"
    start_success = "\033[0;32m"
    start_pass = "\033[0;33m"
    end = "\033[m"
else:
    start_fail = ""
    start_success = ""
    start_pass = ""
    end = ""

#Error code
Exit = 0

