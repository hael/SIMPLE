#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -DGUI:BOOL=off -DCMAKE_BUILD_TYPE:BOOL=debug -DUSE_LIBTIFF:BOOL=on -B build
cmake --build build --config Debug --target clean
