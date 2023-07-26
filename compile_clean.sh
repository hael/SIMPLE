#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DUSE_LIBTIFF=off
make -j install
exit

