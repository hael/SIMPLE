#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DUSE_LIBTIFF=off --warn-uninitialized
make -j install
ctest --parallel 48
#exit

