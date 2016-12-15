
set UNICORE_APP_FOLDER = /opt/unicore/client/UNICORE-Client-7.2.0/applications/
set OPT_APP_FOLDER_LIB = /opt/Application_folder/lib/
set OPT_APP_FOLDER_BIN = /opt/Application_folder/bin/

echo ''
echo '****************************************';
echo 'Now compiling the simple-wrapper...';
echo '****************************************';
echo ''
cd simple-wrapper;
mvn clean install;
cp src/main/resources/bin/simpleWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simpleApplicationWrapper-0.1-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simple-gridbean...';
echo '*****************************************';
echo ''
cd simple-gridbean;
mvn clean install;
cp target/SimpleGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '****************************************';
echo 'Now compiling the simple-preproc-wrapper...';
echo '****************************************';
echo ''
cd simplePreProc-wrapper;
mvn clean install;
cp src/main/resources/bin/simplePreProcWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simplePreProcApplicationWrapper-2.0-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simplepreproc-gridbean...';
echo '*****************************************';
echo ''
cd simplepreproc-gridbean;
mvn clean install;
cp target/SimplePreProcGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '****************************************';
echo 'Now compiling the simplePrime-wrapper...';
echo '****************************************';
echo ''
cd simplePrime-wrapper;
mvn clean install;
cp src/main/resources/bin/simplePrimeWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simplePrimeApplicationWrapper-2.0-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simpleprime-gridbean...';
echo '*****************************************';
echo ''
cd simpleprime-gridbean;
mvn clean install;
cp target/SimplePrimeGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '****************************************';
echo 'Now compiling the simpleEoRecVol-wrapper...';
echo '****************************************';
echo ''
cd simpleEoRecVol-wrapper;
mvn clean install;
cp src/main/resources/bin/simpleEoRecVolWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simpleEoRecVolApplicationWrapper-1.0-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simpleeorecvol-gridbean...';
echo '*****************************************';
echo ''
cd simpleeorecvol-gridbean;
mvn clean install;
cp target/SimpleEoRecVolGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '****************************************';
echo 'Now compiling the simpleEoVolAssemble-wrapper...';
echo '****************************************';
echo ''
cd simpleEoVolAssemble-wrapper;
mvn clean install;
cp src/main/resources/bin/simpleEoVolAssembleWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simpleEoVolAssembleApplicationWrapper-1.0-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simpleeovolassemble-gridbean...';
echo '*****************************************';
echo ''
cd simpleeovolassemble-gridbean;
mvn clean install;
cp target/SimpleEoVolAssembleGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '****************************************';
echo 'Now compiling the simpleStackops-wrapper...';
echo '****************************************';
echo ''
cd simpleStackops-wrapper;
mvn clean install;
cp src/main/resources/bin/simpleStackopsWrapper.sh $OPT_APP_FOLDER_BIN;
cp target/simpleStackopsApplicationWrapper-1.0-SNAPSHOT.jar $OPT_APP_FOLDER_LIB;
cd ../;

echo ''
echo '*****************************************';
echo 'Now compiling the simpleeovolassemble-gridbean...';
echo '*****************************************';
echo ''
cd simplestackops-gridbean;
mvn clean install;
cp target/SimpleStackopsGridBean-2.0.jar $UNICORE_APP_FOLDER;
cd ../;

echo ''
echo '************************************************************';
echo 'end of the compilation of the gridbean and their wrappers...';
echo '************************************************************';
echo ''

