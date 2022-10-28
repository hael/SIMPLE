#!/bin/bash
module purge
module load cmake/3.13.1
echo "Compiling SIMPLE with gcc11.2.0"
module load gcc/11.2.0
rm -r build_gcc11.2.0
mkdir build_gcc11.2.0
cd build_gcc11.2.0
cmake ..
make -j install
cd ..
chmod a+rx,o+rx -R build_gcc11.2.0
module unload gcc/11.2.0
exit

