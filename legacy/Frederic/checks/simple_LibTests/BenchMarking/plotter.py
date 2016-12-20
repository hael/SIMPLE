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
#from mayavi.mlab import *
################################################################################
#                                                                              #
#                                                                              #
#                       Subroutines used in script                             #
#                                                                              #
################################################################################
################################################################################
# Launching the computation                                                    #
################################################################################
def getting_data_from_nrotnk_file(filename):
    print filename
    g_n_b_y_f_n_s_0 = open(filename,'r')
    try:
        g_n_b_y_f_n_s_0_nrot               = range(array_size)
        g_n_b_y_f_n_s_0_nk                 = range(array_size)
        g_n_b_y_f_n_s_0_time               = range(array_size)
        g_n_b_y_f_n_s_0_time_err           = range(array_size)

        iline = -1
        for line in g_n_b_y_f_n_s_0:
            columns = line.split(',')
            columns = [col.strip() for col in columns]
            iline = iline + 1
            from decimal import *
            #print iline,Decimal(columns[0]),Decimal(columns[1]),Decimal(columns[0])/Decimal(columns[1])
            print iline,columns[0], columns[1], columns[2], columns[3]
            g_n_b_y_f_n_s_0_nrot[iline]    = columns[0]
            g_n_b_y_f_n_s_0_nk[iline]      = columns[1]
            g_n_b_y_f_n_s_0_time[iline]    = columns[2]
            g_n_b_y_f_n_s_0_time_err[iline]= columns[3]

            nline = iline + 1
    
    finally:
        g_n_b_y_f_n_s_0.close()

    print "number of lines in gpu_no_bench_yes_fix_no_set_0_z_prime3D_nrotnk101.asc: ",nline
    return
################################################################################
# Launching the computation                                                    #
################################################################################

def graph_xydy_data(io,file_in,xaxis,yaxis,xaxis_lab,yaxis_lab,symbol,line,axes):
    #string constructor
    axes_commands=[xaxis,yaxis]
    axes=' '.join(axes_commands)
    axes_lab_commands=[xaxis_lab,yaxis_lab]
    axes_lab=' '.join(axes_lab_commands)
    graph_commands=[graph,symbol,line,axes,axes_lab,title,io,file_in]
    mygraph=' '.join(graph_commands)
    print mygraph
    os.system('%s'%mygraph)

    return

################################################################################
# Start of the plotting script                                                 #
################################################################################
print "\n"
print "------------------------------------------------------------------------"
print "         Python script to plot output data  the STD and MKL data,       "
print "            By Frederic D.R. Bonnet date: 29th Jul. 2016.               "
print "------------------------------------------------------------------------"
print "\n"

#some parameters settings 
array_size = 250 #maximum size of the arrays
graph = 'graph -T X --bg-color black --frame-color white -I e'
print "Hello python script"
#
#reading the data from the files.
#

############# The STD ##############
#filename = './checks/simple_LibTests/BenchMarking/gpu_no_bench_yes_fix_no_set_0_z_prime3D_nrotnk101.asc'
filename = './ptcls0128/gpu_yes_bench_no_fix_no_set_0_z_prime3D_nrotnk101.asc'
getting_data_from_nrotnk_file(filename)

############# The plotter ##########
#2D plots on the fly using graph
symbol= '--symbol 4'
line='--line-mode 0 -C'
xaxis='-x 0 130'
yaxis='-y 0 2000'
xaxis_lab='--x-label "resolution (nk)"'
yaxis_lab='--y-label "time (secs)"'
io='<'
#getting the files
#GPU
#file_in='./checks/simple_LibTests/BenchMarking/gpu_no_bench_yes_fix_no_set_0_z_prime3D_nk101.asc'
file_in='./ptcls0128/gpu_yes_bench_no_fix_no_set_0_z_prime3D_nk101.asc'
title='-L "Time taken as function of nk GPU-0"'
graph_xydy_data(io,file_in,xaxis,yaxis,xaxis_lab,yaxis_lab,symbol,line,axes)

#file_in='./checks/simple_LibTests/BenchMarking/gpu_yes_bench_no_fix_no_set_0_z_prime3D_nk101.asc'
#title='-L "Time taken as function of nk GPU-2"'
#file_in='./ptcls0128/gpu_yes_bench_no_fix_yes_set_2_z_prime3D_nk101.asc'
#graph_xydy_data(io,file_in,xaxis,yaxis,xaxis_lab,yaxis_lab,symbol,line,axes)

#CPU 
file_in='./ptcls0128/gpu_no_bench_yes_fix_no_set_0_z_prime3D_nk101.asc'
title='-L "Time taken as function of nk bench(CPU)-0"'
graph_xydy_data(io,file_in,xaxis,yaxis,xaxis_lab,yaxis_lab,symbol,line,axes)

file_in='./ptcls0128/gpu_no_bench_no_fix_no_set_0_z_prime3D_nk101.asc'
title='-L "Time taken as function of nk (CPU)-0"'
graph_xydy_data(io,file_in,xaxis,yaxis,xaxis_lab,yaxis_lab,symbol,line,axes)







#3D plots using matplotlib
#plt.figure()
#p1=plt.plot(g_n_b_y_f_n_s_0_nk, g_n_b_y_f_n_s_0_time, 'bo')

#'bo', 'gs', 'ro', 'ys'
#l1 = legend([p1],["8 ptcls bench_gpu=yes"],loc=2)
#gca().add_artist(l1)

#plt.axis([64,3300,0,95.0])
#plt.ylabel('Time (secs)')
#plt.xlabel('n -> [n x n ]')
#plt.title('(Zget_cpu vs. Zsyt_cpu time)')
