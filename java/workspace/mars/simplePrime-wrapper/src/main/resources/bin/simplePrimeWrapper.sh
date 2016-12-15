#!/bin/bash -l

# Wrapper for running simplePrime job on local UNICORE client-server installation
# --------------------------------------------------------------------------
WORKING_DIR=`pwd`
PREFIX=/opt/Application_folder/lib/; export PREFIX

C1=$PREFIX/simplePrimeApplicationWrapper-1.0-SNAPSHOT.jar
C2=$PREFIX/mmm-openmolgrid-common-1.0-SNAPSHOT.jar
C3=$PREFIX/log4j-1.2.14.jar

export CLASSPATH=$C1:$C2:$C3

java -cp $CLASSPATH edu.kit.mmm.wrapper.simplePrime.simplePrimeLauncher \
--workingDirectory=$WORKING_DIR \
--stk=$STK \
--smpd=$SMPD \
--msk=$MSK \
--vol1=$VOL1 \
--vol2=$VOL2 \
--oritab=$ORITAB \
--deftab=$DEFTAB \
--outfile=$OUTFILE \
--kv=$KV \
--cs=$CS \
--fraca=$FRACA \
--frac=$FRAC \
--pgrp=$PGRP \
--ctf=$CTF \
--hp=$HP \
--lp=$LP \
--refine=$REFINE \
--dynlp=$DYNLP \
--noise=$NOISE \
--diversify=$DIVERSIFY \
--eo=$EO \
--norec=$NOREC \
--mw=$MW \
--nstates=$NSTATES \
--startit=$STARTIT \
--trs=$TRS \
--lpstop=$LPSTOP \
--nspace=$NSPACE \
--amsklp=$AMSKLP \
--nnn=$NNN \
--maxits=$MAXITS \
--lpvalid=$LPVALID \
--find=$FIND \
--fstep=$FSTEP \
--fracvalid=$FRACVALID \
--edge=$EDGE \
--dens=$DENS \
--trsstep=$TRSSTEP \
--nvox=$NVOX \
--inner=$INNER \
--time_per_image=$TIME_PER_IMAGE \
--width=$WIDTH \
--nthr=$NTHR


