#!/bin/bash
# Test code (simple_test_exec and production/tests) is skipped unless --compile-tests is given.
# --traps adds floating-point traps and signalling-NaN initialisation (SIMPLE_DEBUG_TRAPS).
BUILD_TESTS=OFF
DEBUG_TRAPS=OFF
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        --traps) DEBUG_TRAPS=ON ;;
        -h|--help) echo "usage: $(basename "$0") [--compile-tests] [--traps]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} -DSIMPLE_DEBUG_TRAPS=${DEBUG_TRAPS} .. -DCMAKE_BUILD_TYPE=debug
make -j || exit $?
# With --compile-tests, the build-time test gate (scripts/run_fast_gate.sh) runs
# between build and install, as in X: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
