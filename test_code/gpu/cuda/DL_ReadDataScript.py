#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Used to extract the data fromt the EM stack for the DeepLearning Proj with
# Group A and B of the Math scholl
# Author: Frederic Bonnet
# Date: 07/09/2016
#----------------------------------------------------------------------------
################################################################################
# Packages imported for the sript                                              #
################################################################################
#System tools
import numpy as np
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
import math
import os
from optparse import OptionParser
#importing own methods
#packages path imports
sys.path.append('./')
#Header package
from header_NW import print_header
from header_NW import print_done
#File package
from file_utils_NW import generate_tar_ball
from file_utils_NW import convert_to_idx3
from file_utils_NW import convert_asc_to_X
#image visualiser package
from PIL import Image
from image_visualiser import visualise
################################################################################
#modify the path to the data here                                              #
################################################################################
path_to_data_good = '/media/frederic/Toshiba_Blue2/Ribo_3262/DeepLearningCase'
path_to_data_badd = '/media/frederic/Toshiba_Blue2/Ribo_3262/DeepLearningCase/Bad_Ptcls_pk_BoxCache'
#modify the filename of the stack here
stk_good1 = 'npart1_nptcls10000_bz240_smpd1.77_msk70_ribo3262.spi'
stk_good2 = 'Good_sumstack_nptcls10975_bz060_smpd708_msk028_negyes_ribo3262.mrc'
#asc_good2 = 'Good_sumstack_nptcls10975_bz060_smpd708_msk028_negyes_ribo3262.asc'
asc_good2 = 'Good_sumstack_nptcls10816_bz060_smpd708_msk028_negyes_ribo3262.asc'
#asc_good2 = 'Good_sumstack_nptcls10816_bz240_smpd177_msk110_negyes_ribo3262.asc'
#sc_good2 = '3_10_Good_sumstack_nptcls10975_bz060_smpd708_msk028_negyes_ribo3262.asc'
stk_badd  = 'Bad_sumstack_nptcls10975_bz060_smpd708_msk028_negyes_ribo3262.mrc'
asc_badd  = 'Bad_sumstack_nptcls10975_bz060_smpd708_msk028_negyes_ribo3262.asc'
################################################################################
# Global variables use through out the script                                  #
################################################################################
#file handling variables
undscr   = "_"           # underscore
sptr   = "/";            # set the name for the separator.
b_sptr = "\\";           # set the name for the separator.
all    = "*";            # unix commands for all files eg: file*
eq     = "=";            # equal sign for file writting
ext    = ".mrc"
local_path = os.getcwd() #get the local path
print_header()           #throw a header
################################################################################
#constructing the file name                                                    #
################################################################################
params_stk_good1=[path_to_data_good , sptr , stk_good1]
params_stk_good2=[path_to_data_good , sptr , stk_good2]
params_asc_good2=[path_to_data_good , sptr , asc_good2]
params_stk_badd=[path_to_data_badd , sptr , stk_badd]
params_asc_badd=[path_to_data_badd , sptr , asc_badd]
GOOD1_STK_IN = ''.join(params_stk_good1)
GOOD2_STK_IN = ''.join(params_stk_good2)
GOOD2_ASC_IN = ''.join(params_asc_good2)
BADD_STK_IN = ''.join(params_stk_badd)
BADD_ASC_IN = ''.join(params_asc_badd)
print "\033[0;91m The local path is        : ", local_path
print "\033[0;92m The GOOD1 stk file name : ", GOOD1_STK_IN
print "\033[0;92m The GOOD2 stk file name : ", GOOD2_STK_IN
print "\033[0;92m The GOOD2 asc file name : ", GOOD2_ASC_IN
print "\033[0;93m The BADD stk file name : ", BADD_STK_IN
print "\033[0;93m The BADD asc file name : ", BADD_ASC_IN
################################################################################
# Writting the data that needs to be written in the correct format such  #
# tensorflow can readin
################################################################################
mx = 60
my = 60
mz = 10816
#data_X = convert_asc_to_X(GOOD2_ASC_IN,dict2['MRC.mx'],dict2['MRC.my'],dict2['MRC.mz'])
data_X = convert_asc_to_X(GOOD2_ASC_IN,mx,my,mz)
#data_X = convert_asc_to_X(GOOD2_ASC_IN,240,240,10816)
#dats_X=convert_asc_to_X(GOOD2_ASC_IN,5,2,3)

#writting and visualising the image file as
file_out = "Good_data_show_raw_nptcls10816_bz060_msk028.png"
visualise(data_X,mx,my,int(np.sqrt(mz)),file_out)

mz = 10975
data_Y = convert_asc_to_X(BADD_ASC_IN,mx,my,mz)
file_out = "Bad_data_show_raw_nptcls10975_bz060_msk028.png"
visualise(data_Y,mx,my,int(np.sqrt(mz)),file_out)

################################################################################
# Now need to write the tensorflow code for the DeepLearning case using data_X #
################################################################################

#TODO: Happy modelling in DeepLearning...

################################################################################
# Printing the finaliser header methodritten in the correct format such        #
################################################################################
print_done()
################################################################################
# End of g the finaliser header methodritten in the correct format such        #
################################################################################
