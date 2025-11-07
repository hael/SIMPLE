#!/bin/bash
export TMPDIR=$HOME/tmp
mkdir $TMPDIR 
chmod u+rwx -R $TMPDIR 
rm -rf build
cmake -G "Unix Makefiles" -D GUI:BOOL=OFF -D NICE=YES -D USE_LIBTIFF:BOOL=ON -B build
cmake --build build --config Release --target clean
