#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DGUI=off -DUSE_LIBTIFF=on -B build
cmake --build build --config Release --target SIMPLE
cpack -G "ZIP"
