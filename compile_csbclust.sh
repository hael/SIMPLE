#!/bin/bash
module purge
module load cmake/3.25.2
echo "Compiling SIMPLE with gcc12.2.0"
module load gcc/12.2.0
rm -r build_gcc12.2.0
mkdir build_gcc12.2.0
cd build_gcc12.2.0
cmake ..
make -j install
cd ..
chmod a+rx,o+rx -R build_gcc12.2.0
module unload gcc/12.2.0
exit

