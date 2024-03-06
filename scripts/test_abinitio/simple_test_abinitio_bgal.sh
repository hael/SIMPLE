#!/bin/bash
if [ $# -eq 0 ]
then
    var_nrestarts=1
else
    var_nrestarts=$1
fi
simple_exec prg=new_project projname=bgal_all &
wait
cd bgal_all
simple_exec prg=import_particles deftab=/mnt/beegfs/elmlund/testbench/bgal/deftab.txt stk=/mnt/beegfs/elmlund/testbench/bgal/sumstack.mrc cs=2.7 fraca=0.1 kv=300 smpd=1.275 &
wait
simple_exec prg=cluster2D ncls=90 mskdiam=180 nparts=4 nthr=10 objfun=euclid projfile=1_import_particles/bgal_all.simple &
wait
simple_exec prg=map_cavgs_selection stk2=2_cluster2D/cavgs_iter016_ranked.mrc ares=70 &
wait
# 'refine=prob'
simple_exec prg=abinitio_3Dmodel pgrp=c1 mskdiam=180 refine=prob nparts=4 nthr=10 projfile=3_selection/bgal_all.simple iterstats=yes nrestarts=$var_nrestarts &
wait