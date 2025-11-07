#!/bin/bash
module purge
module load cmake/3.25.2
echo "Compiling SIMPLE with gcc11.5.0"
module load gcc/11.5.0
rm -r build_gcc11.5.0
cmake -G "Unix Makefiles" -D GUI:BOOL=off -D USE_LIBTIFF:BOOL=on -B build_gcc11.5.0
cmake --build build_gcc11.5.0 --config Release --target clean
chmod -R 777 build_gcc11.5.0
module unload gcc/11.5.0
exit
