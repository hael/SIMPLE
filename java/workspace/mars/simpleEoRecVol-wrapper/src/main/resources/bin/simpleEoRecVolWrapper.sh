#!/bin/bash -l

# Wrapper for running simpleEoRecVol job on local UNICORE client-server installation
# --------------------------------------------------------------------------
WORKING_DIR=`pwd`
PREFIX=/opt/Application_folder/lib/; export PREFIX

C1=$PREFIX/simpleEoRecVolApplicationWrapper-1.0-SNAPSHOT.jar
C2=$PREFIX/mmm-openmolgrid-common-1.0-SNAPSHOT.jar
C3=$PREFIX/log4j-1.2.14.jar

export CLASSPATH=$C1:$C2:$C3

java -cp $CLASSPATH edu.kit.mmm.wrapper.simpleEoRecVol.simpleEoRecVolLauncher \
--workingDirectory=$WORKING_DIR \
--stk=$STK \
--oritab=$ORITAB \
--smpd=$SMPD \
--msk=$MSK \
--kv=$KV \
--frac=$FRAC \
--pgrp=$PGRP \
--mul=$MUL \
--fromp=$FROMP \
--top=$TOP \
--ctf=$CTF \
--fraca=$FRACA \
--cs=$CS \
--deftab=$DEFTAB \
--state=$STATE \
--alpha=$ALPHA \
--mw=$MW \
--amsklp=$AMSKLP \
--edge=$EDGE \
--dens=$DENS \
--box=$BOX \
--inner=$INNER \
--width=$WIDTH \
--nthr=$NTHR \
--part=$PART

