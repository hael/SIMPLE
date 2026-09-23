#!/bin/bash
# Test code (simple_test_exec and production/tests) is skipped unless --compile-tests is given.
BUILD_TESTS=OFF
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) echo "usage: $(basename "$0") [--compile-tests]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
gccversion="15.2.0"
module purge
module load cmake/3.25.2
echo "Compiling SIMPLE with gcc$gccversion"
module load gcc/$gccversion
rm -r build_gcc$gccversion
mkdir build_gcc$gccversion
cd build_gcc$gccversion
cmake -DBUILD_TESTS=${BUILD_TESTS} -D USE_ARCHOPT=OFF ..
make -j || exit $?
# With --compile-tests, the build-time test gate (scripts/run_fast_gate.sh) runs
# between build and install, as in X: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
cd ..
chmod -R 777 build_gcc$gccversion
module unload gcc/$gccversion
exit ${GATE_RC:-0}
