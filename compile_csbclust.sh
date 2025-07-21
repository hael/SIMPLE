#!/bin/bash
module purge
module load cmake/3.25.2
echo "Compiling SIMPLE with gcc11.5.0"
module load gcc/11.5.0
rm -r build_gcc11.5.0
mkdir build_gcc11.5.0
cd build_gcc11.5.0
cmake ..
make -j install
cd ..
chmod -R 777 build_gcc11.5.0
module unload gcc/11.5.0
exit

