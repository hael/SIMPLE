#!/usr/bin/env python

# Edited from EMAN2 / SPARX
# Moves a map according to the matrices output by Chimera as
# a result of a fitting.

# USAGE:  python align_chimera.py volume_to_rotate.mrc 2.0 rotmat.txt
# INPUT1 Volume to rotate
# INPUT2 Sampling distance in Angstroms
# INPUT3 Rotation matrix from Chimera

from EMAN2 import *
from sparx import *
from math import sqrt
import sys

from os import system
from sys import argv
from sys import exit
from string import atof,atoi


def chi2sx(tchi, nx, ny, nz):
        if nx % 2 == 0:
                na = nx/2
        else:
                na = (nx-1)/2
        if ny % 2 == 0:
                nb = ny/2
        else:
                nb = (ny-1)/2
        if nz % 2 == 0:
                nc = nz/2
        else:
                nc = (nz-1)/2
        vcenter = Vec3f(na,nb,nc)
        shift_sx = tchi*vcenter - vcenter
        txlist = tchi.get_matrix()
        txlist[3]  = shift_sx[0]
        txlist[7]  = shift_sx[1]
        txlist[11] = shift_sx[2]
        tsx = Transform(txlist)
        return tsx

def sx2chi(tsx, nx, ny, nz):
        if nx % 2 == 0:
                na = nx/2
        else:
                na = (nx-1)/2
        if ny % 2 == 0:
                nb = ny/2
        else:
                nb = (ny-1)/2
        if nz % 2 == 0:
                nc = nz/2
        else:
                nc = (nz-1)/2
        vcenter = Vec3f(na,nb,nc)
        shift_chi = tsx*(-vcenter) + vcenter
        txlist = tsx.get_matrix()
        txlist[3]  = shift_chi[0]
        txlist[7]  = shift_chi[1]
        txlist[11] = shift_chi[2]
        tchi = Transform(txlist)
        return tchi

def parse_chimera_rotmat( rotmat_file ):
        lines  = open(rotmat_file,'r').readlines()
        rotmat = [0.] * 12
        cnt    = -1
        for line in lines:
                vals = line.strip().split()
                for i in range(len(vals)):
                        cnt += 1
                        rotmat[cnt] = float(vals[i])
        if cnt != 11:
                print 'Not enough values in:', rotmat_file
                exit(-1)
        return rotmat

# Inputs
# map to be fitted
mapf_file   = sys.argv[1]                       # INPUT1 Volume to rotate
mapf        = get_image(mapf_file)
# pixel size of mapf (in Angstrom):
pixf        = float(sys.argv[2])                # INPUT2 Sampling distance in Angstroms
# rotation matrix from chimera
rotmat_file = sys.argv[3]                       # INPUT3 Rotation matrix from Chimera

# size of mapf:
nx = mapf.get_xsize()
ny = mapf.get_ysize()
nz = mapf.get_zsize()

# Rotation matrices
# reference map:
mat_r  = [1.0, 0.0, 0.0, 0.0, \
        0.0, 1.0, 0.0, 0.0, \
        0.0, 0.0, 1.0, 0.0]
# fitted map:
mat_f = parse_chimera_rotmat(rotmat_file)

###########################################################

# change translation units to pixels of map to be fitted:
mat_r[3] /= pixf; mat_f[3] /= pixf
mat_r[7] /= pixf; mat_f[7] /= pixf
mat_r[11] /= pixf; mat_f[11] /= pixf

chi_r = Transform(mat_r)
chi_f = Transform(mat_f)

# relative transformation:
chi_c = chi_r.inverse()*chi_f

sx_c = chi2sx(chi_c,nx,ny,nz)

params_mov = sx_c.get_params("spider")
map_moved  = rot_shift3D(mapf, params_mov["phi"],params_mov["theta"],params_mov["psi"],params_mov["tx"],params_mov["ty"],params_mov["tz"])
print "Rotated volume: transf.spi"
print "PHI  : ",params_mov["phi"]
print "THETA: ",params_mov["theta"]
print "PSI : ",params_mov["psi"]
print "X: ",params_mov["tx"]
print "Y: ",params_mov["ty"]
print "Z: ",params_mov["tz"]
drop_image(map_moved, "transf.spi","s")

# Matrix output
# cmat = sx_c.get_matrix()

# matfile = open("vol-to-ecol_sx.matrix","w")

# # note that shift units are Angstrom:
# for i in range(3):
#         row = "    %9.6f %9.6f %9.6f %8.4f\n" % (cmat[4*i],cmat[4*i+1],cmat[4*i+2],cmat[4*i+3]*pixf)
#         matfile.write(row)

# matfile.close()

exit(0)
