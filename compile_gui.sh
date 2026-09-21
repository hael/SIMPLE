#!/bin/bash
# Test code (simple_test_exec and production/tests) is skipped unless --compile-tests is given.
BUILD_TESTS=OFF
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) echo "usage: $(basename "$0") [--compile-tests]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
export TMPDIR=$HOME/tmp
mkdir $TMPDIR 
chmod u+rwx -R $TMPDIR 
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} .. -D NICE=ON
make -j install
