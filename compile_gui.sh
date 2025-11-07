#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -D GUI:BOOL=on -D USE_LIBTIFF:BOOL=on -B build
cmake --build build --config Release --target clean
