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
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} .. -D USE_COARRAYS=ON
make -j || exit $?
# With --compile-tests, run both the ordinary fast gate and the capability-gated
# two-image coarray smoke before installation. Either failure prevents install.
if [ "$BUILD_TESTS" = ON ]; then
    "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?
    if [ "${GATE_RC:-0}" = 0 ]; then
        echo "-------------------- COARRAY SMOKE TEST --------------------"
        ctest -R '^coarrays$' --no-tests=error --output-on-failure || GATE_RC=$?
    fi
fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
