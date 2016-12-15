#!/bin/bash -l

# Wrapper for running simpleStackops job on local UNICORE client-server installation
# --------------------------------------------------------------------------
WORKING_DIR=`pwd`
PREFIX=/opt/Application_folder/lib/; export PREFIX

C1=$PREFIX/simpleStackopsApplicationWrapper-1.0-SNAPSHOT.jar
C2=$PREFIX/mmm-openmolgrid-common-1.0-SNAPSHOT.jar
C3=$PREFIX/log4j-1.2.14.jar

export CLASSPATH=$C1:$C2:$C3

java -cp $CLASSPATH edu.kit.mmm.wrapper.simpleStackops.simpleStackopsLauncher \
--workingDirectory=$WORKING_DIR \
--stk=$STK \
--stk2=$STK2 \
--oritab=$ORITAB \
--outstk=$OUTSTK \
--deftab=$DEFTAB \
--filetab=$FILETAB \
--smpd=$SMPD \
--kv=$KV \
--cs=$CS \
--fraca=$FRACA \
--nptcls=$NPTCLS \
--frac=$FRAC \
--sym_class=$SYM_CLASS \
--ctf=$CTF \
--ctfsq=$CTFSQ \
--shalgn=$SHALGN \
--roalgn=$ROALGN \
--masscen=$MASSCEN \
--acf=$ACF \
--phrand=$PHRAND \
--vis=$VIS \
--hp=$HP \
--fromp=$FROMP \
--mul=$MUL \
--trs=$TRS \
--lp=$LP \
--snr=$SNR \
--state=$STATE \
--msk=$MSK \
--bin=$BIN \
--top=$TOP \
--nran=$NRAN \
--newbox=$NEWBOX \
--scale=$SCALE \
--hfun=$HFUN \
--norm=$NORM \
--noise_norm=$NOISE_NORM \
--avg=$AVG \
--rankify=$RANKIFY \
--stats=$STATS \
--compare=$COMPARE \
--neg=$NEG \
--merge=$MERGE \
--ft2img=$FT2IMG \
--frameavg=$FRAMEAVG \
--clip=$CLIP \
--box=$BOX \
--inner=$INNER \
--width=$WIDTH \
--outfile=$OUTFILE \
--fraczero=$FRACZERO \
--tres=$TRES \
--split=$SPLIT \
--nthr=$NTHR

