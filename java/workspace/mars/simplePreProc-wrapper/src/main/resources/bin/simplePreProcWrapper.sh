#!/bin/bash -l
# Wrapper for running simplePreProc job on local UNICORE client-server installation
# --------------------------------------------------------------------------
WORKING_DIR=`pwd`
PREFIX=/opt/Application_folder/lib/; export PREFIX

C1=$PREFIX/simplePreProcApplicationWrapper-2.0-SNAPSHOT.jar
C2=$PREFIX/mmm-openmolgrid-common-1.0-SNAPSHOT.jar
C3=$PREFIX/log4j-1.2.14.jar

export CLASSPATH=$C1:$C2:$C3

java -cp $CLASSPATH edu.kit.mmm.wrapper.simplePreProc.simplePreProcLauncher \
--workingDirectory=$WORKING_DIR \
--simple_data_path=$SIMPLE_DATA_PATH \
--nframes=$NFRAMES \
--nconfig=$NCONFIG \
--sph_abe=$SPH_ABE \
--amp_con=$AMP_CON \
--sze_pwr_spc=$SZE_PWR_SPC \
--data_rot=$DATA_ROT \
--target_file=$TARGET_FILE \
--create_dir=$CREATE_DIR \
--file_pre_handler=$FILE_PRE_HANDLER \
--file_organise_dir=$FILE_ORGANISE_DIR \
--unblur=$UNBLUR \
--ctffind=$CTFFIND \
--unblur_dir=$UNBLUR_DIR \
--ctffind_dir=$CTFFIND_DIR \
--file_pst_handler=$FILE_PST_HANDLER \
--have_boxfile=$HAVE_BOXFILE \
--verbose=$VERBOSE \
--delete_dir=$DELETE_DIR

