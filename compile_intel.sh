#!/bin/bash
cmake . -B ./build -DGUI=off -DUSE_LIBTIFF=off -DCMAKE_C_COMPILER=icx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_CXX_COMPILER=icpx
cmake --build 
