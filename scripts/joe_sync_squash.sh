#!/bin/bash
while true; do
 rsync -rltpv /arctica_F3/20161221_172445_pept2_1_9mgml_EC /glusterfs/newstead/PepT2/2016-12-21-Arctica/
Raw
 cd Movies
 find -L ../Raw -iname "*Fractions.mrc" -exec ln -sf {} . \;
 cd ../Micrographs
 for movie in `find -L ../Movies -iname "*Fractions.mrc"`; do
  fileroot=`basename -s ".mrc" $movie | rev | cut -c 11- | rev`
  if [ ! -f $fileroot.txt ]; then
   echo $movie > $fileroot.txt
   sbatch -n 1 -c 20 --nice=1000 -J Unblur --wrap "simple_exec prg=unblur_movies filetab=$fileroot.txt smpd=1.26 nthr=20 fbody=$fileroot exp_time=60 dose_rate=0.44 kv=200 nthr=16 > $fileroot.log"
  fi
 done
 chmod a+rwx *
 cd ..
 sleep 60
done