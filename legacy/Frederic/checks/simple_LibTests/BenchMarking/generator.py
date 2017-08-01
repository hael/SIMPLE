#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Function to generate data for benchMarking using the 1WCM molecule
# Author: Frederic Bonnet
# Date: 03/08/2016
#----------------------------------------------------------------------------
#System tools
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
import numpy as np
from sys import exit
#Imported functions
from input_run_param import get_run_mask, get_molecules_1WCM,get_msk_array_1WCM
from input_run_param import get_boxsz_1WCM, get_diff_smpd_1WCM
from input_run_param import get_diff_boxsz_1WCM
from input_run_param import get_generate_RefVol_1WCM
from input_run_param import get_smpd_val_1WCM
from input_run_param import get_lp_val, get_ptcls_val
################################################################################
# Local variables                                                              #
################################################################################
undscr   = "_"             # underscore
sptr   = "/";           # set the name for the separator.
b_sptr = "\\";          # set the name for the separator.
all    = "*";           # unix commands for all files eg: file*
eq     = "=";           # equal sign for file writting
spc     = " ";           # equal sign for file writting
ext    = ".mrc"
################################################################################
# The generator for the reference volume                                       #
################################################################################
def generate_RefVol(path_data_1WCM,target_path):

    path_1wcm = [path_data_1WCM,"/","1WCM.pdb"]
    PATH_1WCM_IN = ''.join(path_1wcm)

    base_1wcm = "1WCM"
    bz_str = "bz"
    smpd_str = "smpd"
    box = "--box"
    cma = ","
    apix = "--apix"

    for ibox in range(len(get_diff_boxsz_1WCM())):
        print ibox, len(get_diff_boxsz_1WCM())

        out_mrc = [path_data_1WCM,"/",
                   base_1wcm,undscr,bz_str,str(get_diff_boxsz_1WCM()[ibox]),
                   undscr,smpd_str,str(get_diff_smpd_1WCM()[ibox]),ext]
        OUT_MRC_IN = ''.join(out_mrc)

        box_dim = [str(get_diff_boxsz_1WCM()[ibox]),cma,
                   str(get_diff_boxsz_1WCM()[ibox]),cma,
                   str(get_diff_boxsz_1WCM()[ibox])]
        BOX_DIM_IN = ''.join(box_dim)
    
        input_string = [box,BOX_DIM_IN,
                        apix,str(get_diff_smpd_1WCM()[ibox]),
                        PATH_1WCM_IN, OUT_MRC_IN ]

        STRING_IN = ' '.join(input_string)
    
        os.system("ls %s"%PATH_1WCM_IN)

        print "In generate_RefVol STRING_IN: %s"%STRING_IN
        os.system("e2pdb2mrc.py %s"%STRING_IN)

    return
################################################################################
# The generator for the at the command line                                    #
################################################################################
def generate_1WCM(path_data_1WCM,target_path):

    cpu_count = multiprocessing.cpu_count()

    print "def generate_1WCM():"
    print path_data_1WCM
    print target_path
    print cpu_count

    if (get_generate_RefVol_1WCM() in ['yes']):
        generate_RefVol(path_data_1WCM,target_path)
    
    #os.chdir(path_data_1WCM)
    os.system("ls")

    #outer loop for the mask arrays 
    nmsk = -1 
    for imsk in get_msk_array_1WCM(): #[50]:#
        nmsk = nmsk + 1
        #Vol geneeration inner loop for the ipl
        for ilp in get_lp_val(): #['004']:#
            input_string_vol = [path_data_1WCM,str(imsk),
                                str(get_boxsz_1WCM()[nmsk]),
                                str(get_smpd_val_1WCM()[nmsk]),ilp,
                                str(cpu_count)]
            STRING_VOL_IN = ' '.join(input_string_vol)

            print "In generate_1WCM STRING_IN: %s"%STRING_VOL_IN
            os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh generatorData_Vol.csh ./ %s'%STRING_VOL_IN)

        #stack generation inner loop over iptcls    
        for iptcls in get_ptcls_val(): #['0008']:#
            input_string_stk = [path_data_1WCM,str(imsk),
                                str(get_boxsz_1WCM()[nmsk]),
                                str(get_smpd_val_1WCM()[nmsk]),iptcls,
                                str(cpu_count)]
            STRING_STK_IN = ' '.join(input_string_stk)

            print "In generate_1WCM STRING_IN: %s"%STRING_STK_IN
            os.system('cd ./checks/simple_LibTests/BenchMarking/ ; csh generatorData_Stk.csh ./ %s'%STRING_STK_IN)
            
    return
