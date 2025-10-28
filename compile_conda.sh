#!/usr/bin/env bash

#check conda exists
if ! command -v conda >/dev/null 2>&1; then
  echo "conda does not exist."
  exit 1
fi

mkdir -p build

#conda create --prefix `pwd`/build/simple-conda conda-forge::gcc=13.3.0 conda-forge::gxx=13.3.0 conda-forge::gfortran=13.3.0 conda-forge::fftw conda-forge::libtiff conda-forge::jbig conda-forge::jpeg python=3.10 cmake
conda create --prefix `pwd`/build/simple-conda conda-forge::gcc=12.2.0 conda-forge::gxx=12.2.0 conda-forge::gfortran=12.2.0 conda-forge::fftw conda-forge::libtiff conda-forge::jbig conda-forge::jpeg python=3.10 cmake
export PATH=`pwd`/build/simple-conda/bin:$PATH
export LD_LIBRARY_PATH=`pwd`/build/simple-conda/lib:$LD_LIBRARY_PATH
cd build
cmake -D GUI=OFF -D NICE=YES -D TIFF_INCLUDE_DIR=`pwd`/simple-conda/include -D TIFF_LIBRARY_RELEASE=`pwd`/simple-conda/lib/libtiff.so -D CMAKE_PREFIX_PATH=`pwd`/simple-conda ..
make -j install
