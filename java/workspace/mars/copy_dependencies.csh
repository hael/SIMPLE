#!/bin/bash
#ls -R ./*application-wrapper*/src/main/resources/bin/*.sh.plx ./*application-wrapper*/src/main/resources/bin/*.sh >> copy_bin.csh
#/opt/Application_folder/bin

set OPT_APP_FOLDER_LIB = /opt/Application_folder/lib

cp ~/.m2/repository/com/intel/gpe/gpe-client-api/2.2.0/gpe-client-api-2.2.0.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/com/intel/gpe/clients/gpe-gridbean-api/2.2.0/gpe-gridbean-api-2.2.0.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/javaplot/javaplot/0.4.0/javaplot-0.4.0.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/org/jmol/jmol/12.2.3/jmol-12.2.3.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/junit/junit/4.8.2/junit-4.8.2.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/log4j/log4j/1.2.12/log4j-1.2.12.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/log4j/log4j/1.2.14/log4j-1.2.14.jar $OPT_APP_FOLDER_LIB;
cp ~/.m2/repository/org/yaml/snakeyaml/1.10/snakeyaml-1.10.jar $OPT_APP_FOLDER_LIB;


