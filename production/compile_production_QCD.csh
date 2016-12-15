echo "";
echo "*********************************************************";
echo "* Script to compile the production codes                *";
echo "*********************************************************";
echo "";
echo "Now moving back ../ to the ./production directory";
echo "we are in:";
pwd;
######################## Compiling the ./QCD/Impsu3conf/Impsu3conf.f90 ########
echo "Now moving to ./QCD/Impsu3conf directory";
echo "";
cd ./QCD/Impsu3conf;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking Impsu3conf.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;
############ Compiling the ./QCD/Impcoolsu3conf/Impcoolsu3conf_IMPtopo.f90 ####
echo "Now moving to ./QCD/Impcoolsu3conf directory";
echo "";
cd ./QCD/Impcoolsu3conf;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking Impcoolsu3conf_IMPtopo.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;
############ Compiling the ./QCD/Impsmearsu3conf/Impsmearsu3conf_QndImpQ.f90 ##
echo "Now moving to ./QCD/Impsmearsu3conf directory";
echo "";
cd ./QCD/Impsmearsu3conf;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking Impsmearsu3conf_QndImpQ.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;
######################## Compiling the ./QCD/Gaugfixing/Gaugefixsu3conf.f90 ###
echo "Now moving to ./QCD/Gaugfixing/ directory";
echo "";
cd ./QCD/Gaugfixing;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking Gaugefixsu3conf.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;
######################## Compiling the ./QCD/FatlinkClover/FatLinkQP.f90 ######
echo "Now moving to ./QCD/FatlinkClover/ directory";
echo "";
cd ./QCD/FatlinkClover;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking FatLinkQP.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;
##################### Compiling the ./QCD/Massfermion/transform_Jack.f90 ######
##################### Compiling the ./QCD/Massfermion/JackMassAndZ.f90 ########
echo "Now moving to ./QCD/Massfermion/ directory";
echo "";
cd ./QCD/Massfermion;
echo "Cleaning existing compiled/linked code...";
echo "";
make clean;
echo "Compiling and linking transform_Jack.f90 and JackMassAndZ.f90...";
echo "";
make;
echo "Now moving back ../ to the ./production directory";
cd ../../;

##########################

echo "";
echo "*********************************************************";
echo "* The Shell script has been executed                    *";
echo "* The production codes have compiled properly           *";
echo "*********************************************************";
echo "";
echo "Now moving back ../ to the ./production directory";
cd ../../;
