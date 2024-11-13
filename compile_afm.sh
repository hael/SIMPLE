#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DUSE_LIBTIFF=off -DUSE_AFM=ON
make -j install
#exit

