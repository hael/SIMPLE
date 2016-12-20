#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Used to extract the benchmarking results for the SIMPLE library
# Author: Frederic Bonnet
# Date: 23/07/2016
#----------------------------------------------------------------------------
################################################################################
# Packages imported for the sript                                              #
################################################################################
#System tools
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
#Data analysis packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import decimal
import pylab as p
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
#packages path imports
sys.path.append('./checks/simple_LibTests/BenchMarking/')
#Imported functions
from header import print_header
from helper import help
#paths
from input_run_param import get_path_data
from input_run_param import get_path_data_1WCM
from input_run_param import get_target_path
#deciders
from generator import generate_1WCM
from input_run_param import get_simul
from input_run_param import get_parse
from input_run_param import get_generate_1WCM
from input_run_param import get_generate_RefVol_1WCM
#specifiers
from input_run_param import get_maxits
from input_run_param import get_nset
from input_run_param import get_ipart
#Imported functions
from input_run_param import get_lp_val, get_ptcls_val
from input_run_param import get_run_ptcls,get_run_lp,get_run_mask
#trpv1 and pfrib80s
from input_run_param import get_smpd_val_trpv1_pfrid80s, get_boxsz_trpv1_pfrid80s
from input_run_param import get_molecules_trpv1_pfrid80s,get_msk_array_trpv1_pfrid80s
#1WCM
from input_run_param import get_smpd_val_1WCM, get_boxsz_1WCM
from input_run_param import get_molecules_1WCM, get_msk_array_1WCM
#Diff value for the 1WCM data
from input_run_param import get_diff_boxsz_1WCM,get_diff_smpd_1WCM
#full data concanated sets
from input_run_param import get_smpd_val_conc, get_boxsz_conc
from input_run_param import get_molecules_conc, get_msk_array_conc
#chekcer code
from input_run_param import check_conc_arrays
################################################################################
# Global variables use through out the script                                  #
################################################################################
#file handling variables
undscr   = "_"             # underscore
sptr   = "/";           # set the name for the separator.
b_sptr = "\\";          # set the name for the separator.
all    = "*";           # unix commands for all files eg: file*
eq     = "=";           # equal sign for file writting
ext    = ".mrc"
local_path = os.getcwd()
################################################################################
# Global variables use through out the script for file construction and        #
# decision making                                                              #
################################################################################
#file constructors
simul  = get_simul()
parse  = get_parse()
maxits = get_maxits()   #maximum iterations
nset   = get_nset()     #number of experiments 0:gpu={yes,no}, 
ipart  = get_ipart()    #n partition = 1 no distri_pl
#file constructors
npart  = "npart"
volume = "vol"
nptcls = "nptcls"
bz     = "bz"
smpd   = "smpd"
msk    = "msk"
lp     = "lp"
################################################################################
# Local variables used as for the script                                       #
################################################################################
#local variables
path_data      = get_path_data()      #"/media/frederic/Toshiba_Red/BenchMarking"
path_data_1WCM = get_path_data_1WCM() #"/media/frederic/Toshiba_Red/BenchMarking/1WCM"
target_path    = get_target_path()    #"/scratch/frederic/Generated_data"

smpd_val = get_smpd_val_conc()
check_conc_arrays(smpd_val)
#[1.2156, 1.34,
#          1.77, 1.77,
#          1.416, 1.416,
#          1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416,
#          1.011, 1.011, 1.011, 1.011, 1.011, 1.011, 1.011, 1.011, 1.011, 1.011,
#          1.011 ]
boxsz = get_boxsz_conc()
check_conc_arrays(boxsz)
#[256,360,
#       128,128,
#       160,160,
#       256,256,256,256,256,256,256,256,256,256,
#       288,288,288,288,288,288,288,288,288,288,288]

#mol_trpv1,mol_pfrib80s
molecules = get_molecules_conc()
check_conc_arrays(molecules)
#[mol_trpv1, mol_pfrib80s,
#             mol_1WCM128_50, mol_1WCM128_55,
#             mol_1WCM160_65, mol_1WCM160_70,
#             mol_1WCM256_65, mol_1WCM256_70, mol_1WCM256_75,  mol_1WCM256_80,  mol_1WCM256_85,
#             mol_1WCM256_90, mol_1WCM256_95, mol_1WCM256_100, mol_1WCM256_110, mol_1WCM256_115,
#             mol_1WCM288_90,mol_1WCM288_95, mol_1WCM288_100,mol_1WCM288_110,mol_1WCM288_115,
#             mol_1WCM288_120,mol_1WCM288_125,mol_1WCM288_130,mol_1WCM288_135]

msk_array = get_msk_array_conc()
check_conc_arrays(msk_array)
#[70, 128,
#           50,55,
#           65, 70,
#           65, 70, 75, 80, 85, 90, 95, 100, 110, 115,
#           90, 95, 100, 105, 110, 115, 120, 125, 130, 135]
################################################################################
# Run subset of Datasets                                                       #
################################################################################
run_1 = [65, 70, 75, 80, 85]
run_2 = [90, 95, 100, 110]
run_mask = get_run_mask()#[70,128]
run_lp = get_run_lp()#['004','005','006','008','010','012','014','016','018','020']:
run_ptcls = get_run_ptcls()#['0008','0016','0032','0064','0256','0512','1024']:
################################################################################
# Datasets                                                                     #
################################################################################
len_molecules = len(molecules)
#initialising a list for the datasets in question
dataset = []
#constructing the data sets
for imol in range(len(molecules)):
    print molecules[imol]
    path_data_trpv1 = [path_data,sptr,msk,str(msk_array[imol]),undscr,molecules[imol]]
    dataset.insert(imol,''.join(path_data_trpv1))
    print imol,dataset[imol]
################################################################################
# Fixers initialisation                                                        #
################################################################################
USE_GPU_IN = 'no'
BENCH_GPU_IN = 'no'
FIX_GPU_IN = 'no'
SET_GPU_IN = '0'
################################################################################
#                                                                              #
#                                                                              #
#                       Subroutines used in script                             #
#                                                                              #
################################################################################
################################################################################
# Additional methods in header.py and helper.py local files                    #
################################################################################
################################################################################
# Platform details getter                                                      #
################################################################################
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
    return 8#cpu_count
################################################################################
# Getters for the fixers                                                       #
################################################################################
def get_use_gpu():
    if sys.argv[1] in ['--use_gpu=yes']:
        USE_GPU_IN = 'yes'
    if sys.argv[1] in ['--use_gpu=no']:
        USE_GPU_IN = 'no'
    print "in get_use_gpu USE_GPU_IN: ",USE_GPU_IN
    return USE_GPU_IN
def get_bench_gpu():
    if sys.argv[2] in ['--bench_gpu=yes']:
        BENCH_GPU_IN = 'yes'
    if sys.argv[2] in ['--bench_gpu=no']:
        BENCH_GPU_IN = 'no'
    print "in get_bench_gpu BENCH_GPU_IN: ",BENCH_GPU_IN
    return BENCH_GPU_IN
def get_fix_gpu():
    if sys.argv[3] in ['--fix_gpu=yes']:
        FIX_GPU_IN = 'yes'
    if sys.argv[3] in ['--fix_gpu=no']:
        FIX_GPU_IN = 'no'
    print "in get_fix_gpu FIX_GPU_IN: ",FIX_GPU_IN
    return FIX_GPU_IN
def get_set_gpu():
    if sys.argv[4] in ['--set_gpu=0']:
        SET_GPU_IN = '0'
    if sys.argv[4] in ['--set_gpu=1']:
        SET_GPU_IN = '1'
    if sys.argv[4] in ['--set_gpu=2']:
        SET_GPU_IN = '2'
    if sys.argv[4] in ['--set_gpu=3']:
        SET_GPU_IN = '3'
    if sys.argv[4] in ['--set_gpu=4']:
        SET_GPU_IN = '4'
    if sys.argv[4] in ['--set_gpu=5']:
        SET_GPU_IN = '5'
    if sys.argv[4] in ['--set_gpu=6']:
        SET_GPU_IN = '6'
    if sys.argv[4] in ['--set_gpu=7']:
        SET_GPU_IN = '7'
    if sys.argv[4] in ['--set_gpu=8']:
        SET_GPU_IN = '8'
    print "in get_set_gpu SET_GPU_IN: ",SET_GPU_IN
    return SET_GPU_IN
################################################################################
# Launching the computation                                                    #
################################################################################
def launch_calculator(ipart,imol,maxits):

    print "Molecule[imol]: ",molecules[imol]
    print "Dataset[imol]: %s",imol,dataset[imol]
    print "boxsz[imol]: ",imol,dataset[imol]
    print "smpd_val[imol]: ",smpd_val[imol]
    print "maxits: ",maxits

    nchunk = 0
    for iptcls in run_ptcls:
        nchunk = nchunk + 1
        nlp = 0
        for ilp in run_lp:
            nlp = nlp + 1
            nmsk = 0
            for imsk in run_mask:
                nmsk = nmsk + 1
                print "\033[0;93m value of iptcls: ",iptcls," value of ilp",ilp," value of imsk",imsk
                params_stack=[dataset[imol]  ,sptr,
                              npart          ,str(ipart)         ,undscr,
                              nptcls         ,iptcls             ,undscr,
                              bz             ,str(boxsz[imol])   ,undscr,
                              smpd           ,str(smpd_val[imol]),undscr,
                              msk            ,str(imsk)          ,undscr,
                              molecules[imol],ext                      ]
                params_vol=[dataset[imol],sptr,volume,undscr,lp,ilp,undscr,molecules[imol],ext]
                FILE_STACK_IN = ''.join(params_stack)
                FILE_VOL_IN = ''.join(params_vol)
                MOLECULES_IN = molecules[imol]

                ifile = 0 # the condition counter
                
                if os.path.isfile(FILE_STACK_IN) and os.access(FILE_STACK_IN, os.R_OK):
                    print "\033[0;92m File: ",FILE_STACK_IN," exists and is readable"
                    ifile = ifile + 1
                else:
                    print "\033[0;91m Either file: ",FILE_STACK_IN," is missing or is not readable"
                if os.path.isfile(FILE_VOL_IN) and os.access(FILE_VOL_IN, os.R_OK):
                    print "\033[0;92m File: ",FILE_VOL_IN," exists and is readable"
                    ifile = ifile + 1
                else:
                    print "\033[0;91m Either file: ",FILE_VOL_IN," is missing or is not readable"

                print "ifile: ",ifile

                IPTCLS_IN =str(iptcls)
                LP_IN = str(ilp)
                CPU_COUNT_IN = str(cpu_count)
                BZ_IN = str(boxsz[imol])
                SMPD_IN = str(smpd_val[imol])
                IMSK_IN = str(imsk)
                MAXITS_IN = str(maxits)
                input_string = [FILE_STACK_IN,FILE_VOL_IN, MOLECULES_IN,
                                IPTCLS_IN, LP_IN, CPU_COUNT_IN, BZ_IN, SMPD_IN, IMSK_IN,
                                get_use_gpu(), get_bench_gpu(), get_fix_gpu(), get_set_gpu(),MAXITS_IN]
                STRING_IN = ' '.join(input_string)

                input_parse_string = [MOLECULES_IN,
                                      IPTCLS_IN, LP_IN, CPU_COUNT_IN, BZ_IN, SMPD_IN, IMSK_IN,
                                      get_use_gpu(), get_bench_gpu(), get_fix_gpu(), get_set_gpu(),MAXITS_IN]
                STRING_PARSE_IN = ' '.join(input_parse_string)
                
                if (ifile == 2):
                    print "\033[0;92m ******************** %s DataSet... ********************\033[0m"%molecules[imol]
                    #creating and parsing the data
                    if ( simul in ['yes'] and parse in ['yes'] ):
                        os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh create_simulData.csh ./ %s'%STRING_IN)
                        os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh parse_simulData.csh ./ %s'%STRING_PARSE_IN)
                    #creating the data
                    if ( simul in ['yes'] and parse in ['no'] ):
                        os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh create_simulData.csh ./ %s'%STRING_IN)
                    #getting the data from the created data
                    if ( simul in ['no'] and parse in ['yes'] ):
                        os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh parse_simulData.csh ./ %s'%STRING_PARSE_IN)

    #write to file the parameter file generated from the simul
    #spit to file BenchMrkDtAnlzr-incFile.f90
    output = open('./checks/simple_LibTests/BenchMarking/BenchMrkDtAnlzr-incFile.f90','w')
    mxit_nset_string      = ["integer, parameter            ::   nset =",str(nset)]   #n runs with keys use_gpu=yes,no etc
    mxit_nchunk_string    = ["integer, parameter            :: nchunk =",str(nchunk)] #len(['0008']) nspace
    mxit_nlp_string       = ["integer, parameter            ::    nlp =",str(nlp)]    #len(['004','005',..,'020'])
    mxit_nmsk_string      = ["integer, parameter            ::   nmsk =",str(nmsk)]   #len([70,128])
    mxit_MAXITS_IN_string = ["integer, parameter            :: maxits =",str(MAXITS_IN)]
    output.seek(0,2)
    L1 = ' '.join(mxit_nset_string)
    L2 = ' '.join(mxit_nchunk_string)
    L3 = ' '.join(mxit_nlp_string)
    L4 = ' '.join(mxit_nmsk_string)
    L5 = ' '.join(mxit_MAXITS_IN_string)
    output.write(L1+'\n')
    output.write(L2+'\n')
    output.write(L3+'\n')
    output.write(L4+'\n')
    output.write(L5+'\n')
    output.close()
    return
################################################################################
# Launching the analysor                                                       #
################################################################################
def launch_analysor():
    #TODO: here insert the commands for launching the .f90 code
    return
################################################################################
# Launching the analysor                                                       #
################################################################################
def plot_data():
    #import plotter
    #TODO: bring in the 
    return
################################################################################
# Start of benchmarking script                                                 #
################################################################################
print_header()
cpu_count = myplatform()
len_sys = len(sys.argv)
if (len_sys == 2):
    if sys.argv[1] in ['--help']: help()
    exit(0)
if (len_sys > 2):
    if sys.argv[1] in ['--use_gpu=yes'   , '--use_gpu=no']:   print "sys.argv[1]: ",sys.argv[1]
    if sys.argv[2] in ['--bench_gpu=yes' , '--bench_gpu=no']: print "sys.argv[2]: ",sys.argv[2]
    if sys.argv[3] in ['--fix_gpu=yes'   , '--fix_gpu=no']:   print "sys.argv[3]: ",sys.argv[3]

    if sys.argv[4] in ['--set_gpu=0'      , '--set_gpu=0']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=1'      , '--set_gpu=1']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=2'      , '--set_gpu=2']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=3'      , '--set_gpu=3']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=5'      , '--set_gpu=5']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=6'      , '--set_gpu=6']:     print "sys.argv[4]: ",sys.argv[4]
    if sys.argv[4] in ['--set_gpu=7'      , '--set_gpu=7']:     print "sys.argv[4]: ",sys.argv[4]

    if sys.argv[5] in ['--help=yes'      , '--help=no']:      print "sys.argv[5]: ",sys.argv[5]

    if sys.argv[5] in ['--help=yes']: help()

LOCAL_DIR = os.getcwd()
print "\033[0;92m ********************Starting the benchmarking************"

#!!!!!!!!!!!!
# 1WCM DataSet
#!!!!!!!!!!!!
if ( sys.argv[1] not in ['--use_gpu=yes' , '--use_gpu=no', '--bench_gpu=yes' , '--bench_gpu=no', '--fix_gpu=yes' , '--fix_gpu=no', '--set_gpu=yes' , '--set_gpu=no', '--help=yes' , '--help=no'] and sys.argv[2] not in ['--use_gpu=yes' , '--use_gpu=no', '--bench_gpu=yes' , '--bench_gpu=no', '--fix_gpu=yes' , '--fix_gpu=no', '--set_gpu=yes' , '--set_gpu=no', '--help=yes' , '--help=no']):
    print sys.argv[1]
    print sys.argv[2]
    print sys.argv[3]
    print sys.argv[4]
    print sys.argv[5]
    print "\033[0;92m ******************** 1WCM DataSet... ********************\033[0m"
    print "\033[0;92m **no CPU or GPU or BENCH specifier defaulting on CPU ****\033[0m"
    #os.system('cd ./checks/simple_LibTests/1WCM/ ; ln -s spi.polii polii.spi; csh create_simulData.csh ./polii.spi ./ 16 5 %i 240 1.77 76 no no; rm polii.spi'%cpu_count)
else:

    if ( get_generate_1WCM() in ['yes']):
        generate_1WCM(path_data_1WCM,target_path)
    
    for imol in range(len(molecules)):
        launch_calculator(ipart,imol,maxits)

    #launching the data analysor
    launch_analysor()

    #plt data
    plot_data()
    
print get_use_gpu()
print get_bench_gpu()
print get_fix_gpu()
print get_set_gpu()

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

