#!/bin/bash
rm -rf build
mkdir build
cd build
cmake .. -DGUI=off -DUSE_LIBTIFF=off -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_COMPILER=g++
make -j install
#cmake . -B ./build -DGUI=off -DUSE_LIBTIFF=off -DCMAKE_C_COMPILER=icx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_CXX_COMPILER=icpx
#cmake --build 
