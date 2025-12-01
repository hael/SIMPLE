#!/bin/bash
gccversion="15.2.0"
module purge
module load cmake/3.25.2
echo "Compiling SIMPLE with gcc$gccversion"
module load gcc/$gccversion
rm -r build_gcc$gccversion
mkdir build_gcc$gccversion
cd build_gcc$gccversion
cmake ..
make -j install
cd ..
chmod -R 777 build_gcc$gccversion
module unload gcc/$gccversion
exit

