#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off
make -j install
exit

