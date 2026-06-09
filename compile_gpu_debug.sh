#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DUSE_OPENMP_OFFLOAD=ON -DCMAKE_BUILD_TYPE=debug
make -j install
#exit

