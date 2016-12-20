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
from sys import exit
#Data analysis packages
import numpy as np
from getters import *
################################################################################
# Local variables                                                              #
################################################################################
mol_trpv1       = "trpv1"
mol_pfrib80s    = "pfrib80s"

mol_1WCM128_50  = "1WCM128"
mol_1WCM128_55  = "1WCM128"

mol_1WCM160_65  = "1WCM160"
mol_1WCM160_70  = "1WCM160"

mol_1WCM256_65  = "1WCM256"
mol_1WCM256_70  = "1WCM256"
mol_1WCM256_75  = "1WCM256"
mol_1WCM256_80  = "1WCM256"
mol_1WCM256_85  = "1WCM256"
mol_1WCM256_90  = "1WCM256"
mol_1WCM256_95  = "1WCM256"
mol_1WCM256_100 = "1WCM256"
mol_1WCM256_105 = "1WCM256"
mol_1WCM256_110 = "1WCM256"
mol_1WCM256_115 = "1WCM256"

mol_1WCM288_90  = "1WCM288"
mol_1WCM288_95  = "1WCM288"
mol_1WCM288_100 = "1WCM288"
mol_1WCM288_110 = "1WCM288"
mol_1WCM288_115 = "1WCM288"
mol_1WCM288_120 = "1WCM288"
mol_1WCM288_125 = "1WCM288"
mol_1WCM288_130 = "1WCM288"
mol_1WCM288_135 = "1WCM288"
################################################################################
# Run input variables taken into consideration for the                         #
################################################################################
#local variables
#system requirements
#TODO: fix the system requirements
#paths to data
path_data="/media/frederic/Toshiba_Red/BenchMarking/Data" #/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet/checks/simple_LibTests/BenchMarking"
path_data_1WCM="/media/frederic/Toshiba_Red/BenchMarking/Data"#/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet/checks/simple_LibTests/BenchMarking
target_path="/scratch/frederic/Generated_data" #Defined and passed not used yet
#deciders
simul = "no"
parse = "yes"
generate_RefVol_1WCM = "no"
generate_1WCM = "no"
#specifiers
maxits = 6              #maximum iterations
nset   = 2              #number of experiments 0:gpu={yes,no}, 
ipart  = 1              #n partition = 1 no distri_pl
#individual iso runs
run_1    = [65, 70, 75, 80, 85]
run_2    = [90, 95, 100, 110]
run_bench = [50,55,
             65, 70,
             65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115,
             90, 95, 100, 105, 110, 115, 120, 125, 130, 135]
run_cpu   = [50,55,
             65, 70,
             65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115,
             90, 95, 100, 105, 110, 115, 120, 125, 130, 135]
run_set0  = [50,55,
             65, 70,
             65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115,
             90, 95, 100, 105, 110, 115, 120, 125, 130, 135]
#pathches
run_set0_gpu_trpv1 = [50,55]
run_set2_gpu = [130, 135]
run_set2_bench = [95, 100, 105, 110, 115, 120, 125, 130, 135]
run_set3_bench = [125, 130, 135]
run_set2_cpu = [110, 115, 120, 125, 130, 135]
#loop simul runs
run_ptcls= ['0128']#,'0008','0016','0032','0064','0256','0512','1024']
run_lp   = ['004','005','006','008','010','012','014','016','018','020']
#run_mask = run_set0_gpu_trpv1 #run_set3_bench # #run_set0
#run_mask = run_bench
run_mask = run_set0 #run_set0_gpu_trpv1 
################################################################################
# Data range                                                                   #
################################################################################
#The full lp array range
lp_val = ['004', '005', '006', '008', '010', '012', '014', '016', '018', '020']
#The full ptcls array range
ptcls_val = ['0008', '0016', '0032', '0064', '0128', '0256', '0512', '1024']
################################################################################
# trpv1 and pfrib80s Data                                                      #
################################################################################
smpd_val_trpv1_pfrid80s  = [1.2156, 1.34] 
boxsz_trpv1_pfrid80s     = [256,360]
molecules_trpv1_pfrid80s = [mol_trpv1, mol_pfrib80s]
msk_array_trpv1_pfrid80s = [70, 128]
################################################################################
# 1WCM Data                                                                    #
################################################################################
#Reference volume input data
diff_boxsz_1WCM = [128 , 160  , 256  ,   288]
diff_smpd_1WCM =  [1.77, 1.416, 1.416, 1.011]
#full data taken inot consideration
smpd_val_1WCM=[1.77, 1.77,
               1.416, 1.416,
               1.416,1.416,1.416,1.416,1.416,1.416,1.416,1.416,1.416,1.416,1.416,
               1.011,1.011,1.011,1.011,1.011,1.011,1.011,1.011,1.011,1.011,
               1.011]
boxsz_1WCM=[128,128,
            160,160,
            256,256,256,256,256,256,256,256,256,256,256,
            288,288,288,288,288,288,288,288,288,288,288]
molecules_1WCM = [mol_1WCM128_50, mol_1WCM128_55,
                  mol_1WCM160_65, mol_1WCM160_70,
                  mol_1WCM256_65, mol_1WCM256_70, mol_1WCM256_75,
                  mol_1WCM256_80, mol_1WCM256_85, mol_1WCM256_90,
                  mol_1WCM256_95, mol_1WCM256_100, mol_1WCM256_105,
                  mol_1WCM256_110,mol_1WCM256_115,
                  mol_1WCM288_90,mol_1WCM288_95, mol_1WCM288_100,
                  mol_1WCM288_110,mol_1WCM288_115,mol_1WCM288_120,
                  mol_1WCM288_125,mol_1WCM288_130,mol_1WCM288_135]
msk_array_1WCM=[50,55,
                65, 70,
                65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115,
                90, 95, 100, 105, 110, 115, 120, 125, 130, 135]
################################################################################
# The generator for the at the command line (getters)                          #
################################################################################
#calculators
def get_run_mask():
    return run_mask

def get_run_lp():
    return run_lp

def get_run_ptcls():
    return run_ptcls
#full size of the analysed 1WCM data
def get_msk_array_1WCM():
    return msk_array_1WCM

def get_molecules_1WCM():
    return molecules_1WCM

def get_boxsz_1WCM():
    return boxsz_1WCM

def get_smpd_val_1WCM():
    return smpd_val_1WCM

def get_diff_boxsz_1WCM():
    return diff_boxsz_1WCM

def get_diff_smpd_1WCM():
    return diff_smpd_1WCM
#full size of the analysed trpv1 and pfrib80s Data data
def get_smpd_val_trpv1_pfrid80s():
    return smpd_val_trpv1_pfrid80s

def get_boxsz_trpv1_pfrid80s():
    return boxsz_trpv1_pfrid80s

def get_molecules_trpv1_pfrid80s():
    return molecules_trpv1_pfrid80s

def get_msk_array_trpv1_pfrid80s():
    return msk_array_trpv1_pfrid80s
#full size of the analysed trpv1 and pfrib80s Data data
def get_smpd_val_conc():
    return smpd_val_conc

def get_boxsz_conc():
    return boxsz_conc

def get_molecules_conc():
    return molecules_conc

def get_msk_array_conc():
    return msk_array_conc
#paths
def get_path_data():
    return path_data

def get_path_data_1WCM():
    return path_data_1WCM

def get_target_path():
    return target_path

#deciders
def get_generate_1WCM():
    return generate_1WCM

def get_generate_RefVol_1WCM():
    return generate_RefVol_1WCM

def get_simul():
    return simul

def get_parse():
    return parse
#specifiers
def get_maxits():
    return maxits

def get_nset():
    return nset

def get_ipart():
    return ipart
#array rangers
def get_lp_val():
    return lp_val
#The full ptcls array range
def get_ptcls_val():
    return ptcls_val
#checkers code
def check_conc_arrays(conc):
    for i in range(len(conc)):
        print conc[i]
    return
################################################################################
# Combining (by concaneting) to full Data set: 1WCM and trpv1 and pfrib80s     #
################################################################################
#smpd
a = np.array(get_smpd_val_trpv1_pfrid80s())
b = np.array(get_smpd_val_1WCM())
smpd_val_conc = np.concatenate([a, b])
#bxsz
a = np.array(get_boxsz_trpv1_pfrid80s())
b = np.array(get_boxsz_1WCM())
boxsz_conc = np.concatenate([a, b])
#molecules
a = np.array(get_molecules_trpv1_pfrid80s())
b = np.array(get_molecules_1WCM())
molecules_conc = np.concatenate([a, b])
#mask
a = np.array(get_msk_array_trpv1_pfrid80s())
b = np.array(get_msk_array_1WCM())
msk_array_conc = np.concatenate([a, b])
