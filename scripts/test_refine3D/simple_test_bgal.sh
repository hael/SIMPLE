simple_exec prg=new_project projname=bgal &
wait
cd bgal
simple_exec prg=import_particles deftab=/mnt/beegfs/elmlund/testbench/bgal/deftab.txt stk=/mnt/beegfs/elmlund/testbench/bgal/sumstack.mrc cs=2.7 fraca=0.1 kv=300 smpd=1.275 &
wait
simple_exec prg=cluster2D ncls=90 mskdiam=180 nparts=4 nthr=20 objfun=euclid projfile=1_import_particles/bgal.simple &
wait
simple_exec prg=map_cavgs_selection stk2=2_cluster2D/cavgs_iter016_ranked.mrc ares=41 &
wait
simple_exec prg=initial_3Dmodel pgrp=d2 mskdiam=180 nthr=40 projfile=3_selection/bgal.simple &
wait
simple_exec prg=refine3D pgrp=d2 mskdiam=180 nparts=4 nthr=20 maxits=10 refine=neigh objfun=euclid nonuniform=yes sigma_est=global &
wait
simple_exec prg=automask mskdiam=180 amsklp=12 mw=465 thres=0.012 nthr=20 vol1=5_refine3D/recvol_state01_iter010_lp.mrc smpd=1.275 &
wait
simple_exec prg=refine3D pgrp=d2 mskdiam=180 nparts=4 nthr=20 maxits=10 refine=neigh objfun=euclid ml_reg=no nonuniform=yes sigma_est=global continue=yes mskfile=6_automask/automask.mrc combine_eo=yes projfile=5_refine3D/bgal.simple &