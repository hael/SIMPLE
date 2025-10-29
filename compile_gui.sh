#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DNICE=YES -DUSE_LIBTIFF=on
make -j install
#exit

