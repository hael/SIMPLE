#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DCMAKE_BUILD_TYPE=debug
make -j install
exit

