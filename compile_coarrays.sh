#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -D USE_COARRAYS=ON
make -j install
#exit

