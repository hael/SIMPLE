#!/bin/sh

if [ -z "${SIMPLE_TESTBENCH_DATA+set}" ]; then
   echo "SIMPLE_TESTBENCH_DATA env variable is not set."; exit 1
fi

if [ ! -d "${SIMPLE_TESTBENCH_DATA}" ]; then
    echo "SIMPLE_TESTBENCH_DATA variable is not a directory."; exit 1
fi

if [ ! -d "${SIMPLE_TESTBENCH_DATA}/stars_from_matt" ]; then
    echo "stars_from_matt is not a directory in the benchmark data folder."; exit 1
fi
SECONDS=0

set +ev

echo " Creating new project in SimpleStarImport"
[ -d SimpleStarImport ] && rm -rf SimpleStarImport
simple_exec prg=new_project projname=SimpleStarImport
if [ ! -d SimpleStarImport ];then echo new_project failed; exit 1;fi


cd SimpleStarImport

echo " Testing import_starproject without startype and smpd arguments"
simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star
if [ $? -ne 0 ];then
    echo  'EXPECTED:   prg=importstar should fail without startype and smpd';
else
    echo  'UNEXPECTED: prg=importstar should fail without startype and smpd';
    exit 1;
fi

echo " Testing import_starproject without startype argument"
simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star smpd=1.1
if [ $? -ne 0 ];then
    echo 'EXPECTED: prg=importstar should fail without startype ';
else
    echo 'UNEXPECTED: prg=importstar should fail without startype and smpd';
    exit 1;
fi

echo " Testing import_starproject with irregular startype argument"
simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star smpd=1.1 startype=blah
if [ $? -ne 0 ];then
    echo 'EXPECTED: prg=importstar should fail without a valid startype ';
else
    echo 'UNEXPECTED: prg=importstar should fail without a valid startype';
    exit 1;
fi

echo " Testing import_starproject with old-format startype argument"
simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star smpd=1.1 startype=extract
if [ $? -ne 0 ];then
    echo 'EXPECTED: prg=importstar should fail without a valid startype (extract is an older version) ';
else
    echo 'UNEXPECTED: prg=importstar should fail without a valid startype';
    exit 1;
fi

echo " Testing import_starproject with old-format startype argument"
simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star smpd=1.1 startype=extract oritab=oritab-stardoc.txt
if [ $? -ne 0 ];then
    echo 'EXPECTED: prg=importstar should fail ';
else
    echo 'UNEXPECTED: prg=importstar should fail without a valid startype';
    exit 1;
fi


simple_exec prg=import_starproject starfile=${SIMPLE_TESTBENCH_DATA}/stars_from_matt/Extract/364Box_Extract_LocalCTF/particles.star smpd=1.1 startype=particles
if [ $? -eq 0 ];then
    echo 'EXPECTED  valid import star args ';
else
    echo 'UNEXPECTED prg=importstar should not fail with a valid startype and smpd';
    exit 1;
fi

echo '$0 completed successfully '

exit 0;
