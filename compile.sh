#!/bin/bash
rm -rf build
cmake -G "Unix Makefiles" -B build -DGUI=off -DUSE_LIBTIFF=on
cmake --build build --config Release 
