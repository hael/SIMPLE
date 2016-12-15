#compiling the openmolgrid package
cd mmm-openmolgrid-client;
mvn clean install;
cd ../;
cd common;
mvn clean install;
cd ../;
#compiling the gridbeans
cd simple-gridbean;
mvn clean install;
cd ../;
cd simple-wrapper;
mvn clean install;
cd ../;
cd simplepreproc-gridbean;
mvn clean install;
cd ../;
cd simplePreProc-wrapper;
mvn clean install;
cd ../;
cd simpleprime-gridbean;
mvn clean install;
cd ../;
cd simplePrime-wrapper;
mvn clean install;
cd ../;
cd simpleeorecvol-gridbean;
mvn clean install;
cd ../;
cd simpleEoRecVol-wrapper;
mvn clean install;
cd ../;
cd simpleeovolassemble-gridbean;
mvn clean install;
cd ../;
cd simpleEoVolAssemble-wrapper;
mvn clean install;
cd ../;
cd simplestackops-gridbean;
mvn clean install;
cd ../;
cd simpleStackops-wrapper;
mvn clean install;
cd ../;





