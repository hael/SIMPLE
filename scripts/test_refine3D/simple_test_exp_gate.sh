simple_exec prg=new_project projname=exp_gate &
wait
cd exp_gate
simple_exec prg=import_particles cs=2.7 fraca=0.1 kv=300 smpd=1.72 deftab=/mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart.txt stk=/mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart_sc.mrc &
wait
simple_exec prg=cluster2D ncls=25 objfun=euclid sigma_est=global mskdiam=155 nthr=20 nparts=4 projfile=1_import_particles/exp_gate.simple kweight=all &
wait
simple_exec prg=initial_3Dmodel pgrp=c1 mskdiam=155 nthr=40 projfile=2_cluster2D/exp_gate.simple &
wait
simple_exec prg=refine3D pgrp=c1 mskdiam=155 nparts=4 nthr=20 maxits=10 objfun=euclid refine=neigh nonuniform=yes sigma_est=global projfile=3_initial_3Dmodel/exp_gate.simple &
wait
simple_exec prg=automask mskdiam=155 amsklp=15 mw=300 thres=0.055 nthr=20 vol1=4_refine3D/recvol_state01_iter010_lp.mrc smpd=1.72 &
wait
simple_exec prg=refine3D pgrp=c1 mskdiam=155 nparts=4 nthr=20 maxits=10 refine=neigh objfun=euclid nonuniform=yes sigma_est=global continue=yes mskfile=5_automask/automask.mrc projfile=4_refine3D/exp_gate.simple &