#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=on -DUSE_LIBTIFF=on
make -j install
#exit

