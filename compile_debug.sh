#!/bin/bash
# Test programs (production/tests) are skipped unless --compile-tests is given.
BUILD_TESTS=OFF
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) echo "usage: $(basename "$0") [--compile-tests]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} .. -DCMAKE_BUILD_TYPE=debug
make -j install
#exit

