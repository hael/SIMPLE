#!/usr/bin/env bash
# Tests are built by default (simple_test_exec and the *_tester modules) and the fast
# gate runs before installation; --exclude-tests builds the library and executables only.
BUILD_TESTS=ON
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        -h|--help) echo "usage: $(basename "$0") [--exclude-tests]"; exit 0 ;;
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
make -j || exit $?
# Unless --exclude-tests is given, the build-time test gate (scripts/run_fast_gate.sh)
# runs between build and install, as in X: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
