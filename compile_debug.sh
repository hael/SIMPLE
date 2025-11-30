#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DCMAKE_BUILD_TYPE=debug -DUSE_LIBTIFF=on
make -j install
#exit

