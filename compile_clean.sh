#!/bin/bash
# Clean release build and install into build/.
#
#   ./compile_clean.sh                  library and executables only (default)
#   ./compile_clean.sh --compile-tests  also build the test programs in production/tests
#
# The test programs add ~110 s of compile CPU and one static link each; they
# are rarely needed for everyday work, so they are off unless asked for.
BUILD_TESTS=OFF
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "compile_clean.sh: unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake .. -DBUILD_TESTS=${BUILD_TESTS}
make -j install
#exit
