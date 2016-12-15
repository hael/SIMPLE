#!/bin/bash -l

# Script for running DPoly job on CINECA PLX
# ------------------------------------------

C1=/plx/userexternal/sbozic00/m3hpc/mmm_wrapper/lib/simpleStackopsApplicationWrapper-0.1-SNAPSHOT.jar
C2=/plx/userexternal/sbozic00/m3hpc/mmm_wrapper/lib/mmm-openmolgrid-common-0.5-SNAPSHOT.jar
C3=/plx/userexternal/sbozic00/m3hpc/mmm_wrapper/lib/log4j-1.2.14.jar

export CLASSPATH=$C1:$C2:$C3

module load autoload dl_poly

java -Xmx1024m -cp $CLASSPATH edu.kit.mmm.wrapper.simpleStackops.simpleStackopsLauncher \
--workingDirectory=${PWD} \
--mpiCommand="mpirun" \
--numberProcessors=$UC_PROCESSORS
