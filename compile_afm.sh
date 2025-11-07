#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -D USE_LIBTIFF:BOOL=ON -D USE_LIBTIFF:BOOL=OFF -D USE_AFM:BOOL=ON -B build
cmake --build build --config Release --target clean
