#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -D NICE=YES -D USE_LIBTIFF=ON -D GUI=OFF
make -j install
