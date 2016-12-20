#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Helper function to for file utils for the data tranformer for
# DL_ReadDataScript.py
# Author: Frederic Bonnet
# Date: 07/09/2016
#----------------------------------------------------------------------------
#System tools
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
#computation stuff
import numpy as np
from decimal import *
#
################################################################################
# The header subroutine for the script                                         #
################################################################################
import tarfile
################################################################################
# parser from the asc files                                                    #
################################################################################
#simple_stackops stk=flipped.mrc top=XXX outstk=XXX_etc.mrc
def convert_asc_to_X(file_in,nx,ny,nptcls):
    print "\033[0;92m  Box size nx:",nx
    print "\033[0;92m  Box size ny:",ny
    print "\033[0;92m There are nz:",nptcls," in method convert_asc_to_X"
    print "\033[0;93m Looking at file: ",file_in
    file_X = open(file_in,'r')
    lines_X = file_X.readlines()
    file_X.close()
    filename_out = "out.asc"
    #file_out = open(filename_out,'w')
    try:

        #nptcls = 4
        print "\033[0;92m The value of nptcls: ",nptcls

        line = nptcls
        rows = nx*ny
        data_X = np.zeros(shape=(line,rows))
        print "\033[0;93m data_X after initialisation and being populated"
        print data_X
        print "\033[0;92m Number of particles      :",line
        print "\033[0;92m length of rows           :",rows
        print "\033[0;95m Populating data_X... from: ",file_in.strip()
        iline = -1
        for ilines_X in lines_X:
            iline += 1
            vals   = ilines_X.strip().split()
            #print ilin, len(vals), vals
            for irow in range(rows):
                data_X[iline,irow] = float(vals[irow])
                #file_out.write(str(data_X[iline][irow])+' ')
            #file_out.write('\n')

    finally:
        file_X.close()
        #file_out.close()
        print "The 2D matrix data_X has been populated from: ",file_in
    return data_X
################################################################################
# tar ball generator (may need adaptation)                                     #
################################################################################
def generate_tar_ball(tar_file_out, files_to_tar):
    tar = tarfile.open(tar_file_out, "w:gz")
    for name in files_to_tar:    #example ["foo", "bar", "quux"]
        tar.add(name)
        tar.close()
    return
################################################################################
# IDX3 convertor                                                               #
################################################################################

#TODO: here write the method to convert the files into idx3 format

def convert_to_idx3():
    #TODO: insert code here
    return
