#!/bin/bash
export TMPDIR=$HOME/tmp
mkdir $TMPDIR 
chmod u+rwx -R $TMPDIR 
rm -rf build
mkdir build
cd build
cmake .. -D NICE=ON
make -j install
