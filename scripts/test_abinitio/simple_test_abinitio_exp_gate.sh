#!/bin/bash
if [ $# -eq 0 ]
then
    var_nrestarts=1
else
    var_nrestarts=$1
fi

# 'refine=prob'
simple_exec prg=new_project projname=exp_gate_prob &
wait
cd exp_gate_prob
wait
simple_exec prg=stackops fromp=1 top=20000 stk=/mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart_sc.mrc nthr=20 smpd=1.72 &
wait
head -n 20000 /mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart.txt > new_deftab.txt &
wait
simple_exec prg=import_particles cs=2.7 fraca=0.1 kv=300 smpd=1.72 deftab=new_deftab.txt stk=outstk.mrc &
wait
simple_exec prg=abinitio_3Dmodel mskdiam=155 pgrp=c1 refine=prob nparts=4 nthr=10 projfile=1_import_particles/exp_gate_prob.simple iterstats=yes nrestarts=$var_nrestarts &
wait
cd ..
wait

# 'refine=shc_smpl'
simple_exec prg=new_project projname=exp_gate_shc_smpl &
wait
cd exp_gate_shc_smpl
wait
simple_exec prg=stackops fromp=1 top=20000 stk=/mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart_sc.mrc nthr=20 smpd=1.72 &
wait
head -n 20000 /mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart.txt > new_deftab.txt &
wait
simple_exec prg=import_particles cs=2.7 fraca=0.1 kv=300 smpd=1.72 deftab=new_deftab.txt stk=outstk.mrc &
wait
simple_exec prg=abinitio_3Dmodel mskdiam=155 pgrp=c1 refine=shc_smpl nparts=4 nthr=10 projfile=1_import_particles/exp_gate_shc_smpl.simple iterstats=yes nrestarts=$var_nrestarts &
wait
cd ..
wait

# 'refine=prob prob_sh=yes'
simple_exec prg=new_project projname=exp_gate_prob_sh &
wait
cd exp_gate_prob_sh
wait
simple_exec prg=stackops fromp=1 top=20000 stk=/mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart_sc.mrc nthr=20 smpd=1.72 &
wait
head -n 20000 /mnt/beegfs/elmlund/testbench/exp_gate/Refine244_centredpart.txt > new_deftab.txt &
wait
simple_exec prg=import_particles cs=2.7 fraca=0.1 kv=300 smpd=1.72 deftab=new_deftab.txt stk=outstk.mrc &
wait
simple_exec prg=abinitio_3Dmodel mskdiam=155 pgrp=c1 refine=prob nparts=4 nthr=10 projfile=1_import_particles/exp_gate_prob_sh.simple iterstats=yes prob_sh=yes nrestarts=$var_nrestarts &
wait
cd ..
wait