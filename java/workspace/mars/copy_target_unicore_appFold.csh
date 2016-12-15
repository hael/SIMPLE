#!/bin/bash

set APP_FLDR_LIB = /opt/Application_folder/lib
set UNICR_CLT_APP_FLDR = /opt/unicore/client/UNICORE-Client-7.2.0/applications

echo "The path for the APP_FLDR_LIB: $APP_FLDR_LIB"
echo "The path for the UNICR_CLT_APP_FLDR: $UNICR_CLT_APP_FLDR "

#$APP_FLDR_LIB
#$UNICR_CLT_APP_FLDR

################################
#copying the openmolgrid package
################################

#./mmm-openmolgrid-client/target:

#cp ./mmm-openmolgrid-client/target/mmmm-openmolgrid-client-1.0-SNAPSHOT.jar

#./common/target:

cp ./common/target/mmm-openmolgrid-common-1.0-SNAPSHOT.jar $APP_FLDR_LIB;

################################
#copying the GridBeans and the application wrappers
################################

#./simple-gridbean

cp ./simple-gridbean/target/SimpleGridBean-2.0.jar $UNICR_CLT_APP_FLDR;

#./simple-wrapper

cp ./simple-wrapper/target/simpleApplicationWrapper-0.1-SNAPSHOT.jar $APP_FLDR_LIB;

#./simple-preproc-gridbean

cp ./simple-preproc-gridbean/target/SimplePreProcGridBean-2.0.jar $UNICR_CLT_APP_FLDR;

#./simple-preproc-wrapper

cp ./simple-preproc-wrapper/target/simplePreProcApplicationWrapper-2.0-SNAPSHOT.jar $APP_FLDR_LIB;


