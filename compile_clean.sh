#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -D CMAKE_BUILD_TYPE:BOOL=Release -D GUI:BOOL=off -D USE_LIBTIFF:BOOL=on -B build
cmake --build build --config Release --target SIMPLE
cmake --install ./build --prefix /path/to/somewhere
