#!/bin/bash

set MARS = "/home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet/java/workspace/mars"

set changer = "chmod a-x Makefile_target *.f *.f90 Makefile_* old_Makefile_target  *.cpp *.cu *.cuh *.cc *.gz *.txt  *.py *.c *.pdf *.tex *.png *.ai *.docx *.h *macros* junk_* README *.m4 *.in LICENSE Makefile.* config.guess  config.sub  depcomp  *-sh  ltmain.sh  missing *.cfg html/*.html *.html *.css *.dsp  *.dsw *.vcproj *.sln  Makefile  *.m  *.asc *.eps *.jpg *.agr *.ps *.xml *.java *.sh *.class *.awk *.xls *.xyz *.cml *.pdb *.mop *.mopout pdt_15_R.out *.out *.yaml *.ref *.gif *.form *.xsl heme opi1 water *.kpt *.local *.louhi *.plx *.properties simpleidb"

#The changer method
cd ./; $changer;
cd $MARS;

cd ./simple-wrapper; $changer;
cd $MARS;

cd ./simple-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simple-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./common; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/bigdft; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/cml; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/csv; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/gamess; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/mmtk; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/mopac; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/pdb; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/sdf; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/slf; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/smiles; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/tofet; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/xml; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/xyz; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/format/yaml; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/model; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/process; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/qsar; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/util; $changer;
cd $MARS;

cd ./common/src/main/java/org/openmolgrid/wrapper; $changer;
cd $MARS;

cd ./common/src/main/resources/atoms; $changer;
cd $MARS;

cd ./common/src/main/resources/dictionary; $changer;
cd $MARS;

cd ./common/src/test/java/org/openmolgrid/common; $changer;
cd $MARS;

cd ./common/src/test/java/org/openmolgrid/util; $changer;
cd $MARS;

cd ./common/src/test/resources/bigDFT; $changer;
cd $MARS;

cd ./common/src/test/resources/cml; $changer;
cd $MARS;

cd ./common/src/test/resources/deposit; $changer;
cd $MARS;

cd ./common/src/test/resources/mopac; $changer;
cd $MARS;

cd ./common/src/test/resources/xyz; $changer;
cd $MARS;

cd ./common/src/test/resources/yaml; $changer;
cd $MARS;

cd ./common/src/test/resources/yaml/Davidson-SiH4; $changer;
cd $MARS;

cd ./common/src/test/resources/yaml/Li+; $changer;
cd $MARS;

cd ./common/src/test/resources/yaml/O2-Spin; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/main/java/org/chemomentum/gridbeans/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/main/java/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/main/java/org/openmolgrid/client/common/table; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/main/resources/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/test/java/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/bin/src/test/resources/yaml; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/main/java/org/chemomentum/gridbeans/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/main/java/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/main/java/org/openmolgrid/client/common/table; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/main/resources/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/test/java/org/openmolgrid/client/common; $changer;
cd $MARS;

cd ./mmm-openmolgrid-client/src/test/resources/yaml; $changer;
cd $MARS;

cd ./simpleeorecvol-gridbean; $changer;
cd $MARS;

cd ./simpleeorecvol-gridbean/src/main/java/edu/kit/gridbeans/simpleeorecvol; $changer;
cd $MARS;

cd ./simpleeorecvol-gridbean/src/main/java/edu/kit/gridbeans/simpleeorecvol/plugin; $changer;
cd $MARS;

cd ./simpleeorecvol-gridbean/src/main/resources/edu/kit/gridbeans/simpleeorecvol/images; $changer;
cd $MARS;

cd ./simpleeorecvol-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simpleEoRecVol-wrapper; $changer;
cd $MARS;

cd ./simpleEoRecVol-wrapper/src/main/java/edu/kit/mmm/wrapper/simpleEoRecVol; $changer;
cd $MARS;

cd ./simpleEoRecVol-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simpleEoRecVol-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simpleeovolassemble-gridbean; $changer;
cd $MARS;

cd ./simpleeovolassemble-gridbean/src/main/java/edu/kit/gridbeans/simpleeovolassemble; $changer;
cd $MARS;

cd ./simpleeovolassemble-gridbean/src/main/java/edu/kit/gridbeans/simpleeovolassemble/plugin; $changer;
cd $MARS;

cd ./simpleeovolassemble-gridbean/src/main/resources/edu/kit/gridbeans/simpleeovolassemble/images; $changer;
cd $MARS;

cd ./simpleeovolassemble-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simpleEoVolAssemble-wrapper; $changer;
cd $MARS;

cd ./simpleEoVolAssemble-wrapper/src/main/java/edu/kit/mmm/wrapper/simpleEoVolAssemble; $changer;
cd $MARS;

cd ./simpleEoVolAssemble-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simpleEoVolAssemble-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simple-gridbean; $changer;
cd $MARS;

cd ./simple-gridbean/src/main/java/edu/kit/gridbeans/simple; $changer;
cd $MARS;

cd ./simple-gridbean/src/main/java/edu/kit/gridbeans/simple/plugin; $changer;
cd $MARS;

cd ./simple-gridbean/src/main/resources/edu/kit/gridbeans/simple/images; $changer;
cd $MARS;

cd ./simple-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simplepreproc-gridbean; $changer;
cd $MARS;

cd ./simplepreproc-gridbean/src/main/java/edu/kit/gridbeans/simplepreproc; $changer;
cd $MARS;

cd ./simplepreproc-gridbean/src/main/java/edu/kit/gridbeans/simplepreproc/plugin; $changer;
cd $MARS;

cd ./simplepreproc-gridbean/src/main/resources/edu/kit/gridbeans/simplepreproc/images; $changer;
cd $MARS;

cd ./simplepreproc-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simplePreProc-wrapper; $changer;
cd $MARS;

cd ./simplePreProc-wrapper/src/main/java/edu/kit/mmm/wrapper/simplePreProc; $changer;
cd $MARS;

cd ./simplePreProc-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simplePreProc-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simpleprime-gridbean; $changer;
cd $MARS;

cd ./simpleprime-gridbean/src/main/java/edu/kit/gridbeans/simpleprime; $changer;
cd $MARS;

cd ./simpleprime-gridbean/src/main/java/edu/kit/gridbeans/simpleprime/plugin; $changer;
cd $MARS;

cd ./simpleprime-gridbean/src/main/resources/edu/kit/gridbeans/simpleprime/images; $changer;
cd $MARS;

cd ./simpleprime-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simplePrime-wrapper; $changer;
cd $MARS;

cd ./simplePrime-wrapper/src/main/java/edu/kit/mmm/wrapper/simplePrime; $changer;
cd $MARS;

cd ./simplePrime-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simplePrime-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simplestackops-gridbean; $changer;
cd $MARS;

cd ./simplestackops-gridbean/src/main/java/edu/kit/gridbeans/simplestackops; $changer;
cd $MARS;

cd ./simplestackops-gridbean/src/main/java/edu/kit/gridbeans/simplestackops/plugin; $changer;
cd $MARS;

cd ./simplestackops-gridbean/src/main/resources/edu/kit/gridbeans/simplestackops/images; $changer;
cd $MARS;

cd ./simplestackops-gridbean/src/main/resources/META-INF; $changer;
cd $MARS;

cd ./simpleStackops-wrapper; $changer;
cd $MARS;

cd ./simpleStackops-wrapper/src/main/java/edu/kit/mmm/wrapper/simpleStackops; $changer;
cd $MARS;

cd ./simpleStackops-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simpleStackops-wrapper/src/main/resources/bin; $changer;
cd $MARS;

cd ./simple-wrapper; $changer;
cd $MARS;

cd ./simple-wrapper/src/main/java/edu/kit/mmm/wrapper/simple; $changer;
cd $MARS;

cd ./simple-wrapper/src/main/resources; $changer;
cd $MARS;

cd ./simple-wrapper/src/main/resources/bin; $changer;
cd $MARS;
