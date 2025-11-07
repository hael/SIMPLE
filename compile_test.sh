#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -DGUI:BOOL=off -DUSE_LIBTIFF:BOOL=off --warn-uninitialized -B build
cmake --build build --config Release --target clean
make -j install
ctest --parallel 48
