#!/bin/bash
#ls -R ./*application-wrapper*/src/main/resources/bin/*.sh.plx ./*application-wrapper*/src/main/resources/bin/*.sh >> copy_bin.csh
#/opt/Application_folder/bin

set OPT_APP_FOLDER_BIN = /opt/Application_folder/bin

cp ./simple-wrapper/src/main/resources/bin/simpleWrapper.sh $OPT_APP_FOLDER_BIN;
cp ./simplePreProc-wrapper/src/main/resources/bin/simplePreProcWrapper.sh $OPT_APP_FOLDER_BIN;
cp ./simpleEoRecVol-wrapper/src/main/resources/bin/simpleEoRecVolWrapper.sh $OPT_APP_FOLDER_BIN;
cp ./simpleEoVolAssemble-wrapper/src/main/resources/bin/simpleEoVolAssembleWrapper.sh $OPT_APP_FOLDER_BIN;
cp ./simpleStackops-wrapper/src/main/resources/bin/simpleStackopsWrapper.sh $OPT_APP_FOLDER_BIN;
cp ./simplePrime-wrapper/src/main/resources/bin/simplePrimeWrapper.sh $OPT_APP_FOLDER_BIN;
