#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -D USE_COARRAYS=on
make -j install
#exit

