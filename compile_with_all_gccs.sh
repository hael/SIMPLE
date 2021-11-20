#!/bin/bash
module purge
module load cmake/3.13.1

echo "Compiling SIMPLE with gcc5.5.0"
module load gcc/5.5.0
rm -r build_gcc5.5.0
mkdir build_gcc5.5.0
cd build_gcc5.5.0
cmake ..
make -j 8 install
cd ..
module unload gcc/5.5.0

echo "Compiling SIMPLE with gcc9.3.0"
module load gcc/9.3.0
rm -r build_gcc9.3.0
mkdir build_gcc9.3.0
cd build_gcc9.3.0
cmake ..
make -j 8 install
cd ..
module unload gcc/9.3.0

echo "Compiling SIMPLE with gcc10.2.0"
module load gcc/10.2.0
rm -r build_gcc10.2.0
mkdir build_gcc10.2.0
cd build_gcc10.2.0
cmake ..
make -j 8 install
cd ..
module unload gcc/10.2.0

exit
