#!/usr/bin/env bash
# Test programs (production/tests) are skipped unless --compile-tests is given.
BUILD_TESTS=OFF
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) echo "usage: $(basename "$0") [--compile-tests]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done

#check conda exists
if ! command -v conda >/dev/null 2>&1; then
  echo "conda does not exist."
  exit 1
fi

mkdir -p build

conda create --prefix `pwd`/build/simple-conda conda-forge::gcc=15.2.0 conda-forge::gxx=15.2.0 conda-forge::gfortran=15.2.0 conda-forge::fftw conda-forge::libtiff conda-forge::jbig conda-forge::jpeg python=3.10 cmake
export PATH=`pwd`/build/simple-conda/bin:$PATH
export LD_LIBRARY_PATH=`pwd`/build/simple-conda/lib:$LD_LIBRARY_PATH
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} -D NICE=ON -D TIFF_INCLUDE_DIR=`pwd`/simple-conda/include -D TIFF_LIBRARY_RELEASE=`pwd`/simple-conda/lib/libtiff.so -D CMAKE_PREFIX_PATH=`pwd`/simple-conda ..
make -j install
