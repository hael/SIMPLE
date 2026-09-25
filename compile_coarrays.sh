#!/bin/bash
# Tests are built by default (simple_test_exec and the *_tester modules) and the fast
# gate runs before installation; --exclude-tests builds the library and executables only.
# High-level tests are registered but run only by an explicit `ctest -L highlevel`.
BUILD_TESTS=ON
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        -h|--help) echo "usage: $(basename "$0") [--exclude-tests]"; exit 0 ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} .. -D USE_COARRAYS=ON
make -j || exit $?
# Unless --exclude-tests is given, run both the ordinary fast gate and the
# capability-gated two-image coarray smoke before installation. Either failure prevents install.
if [ "$BUILD_TESTS" = ON ]; then
    "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?
    if [ "${GATE_RC:-0}" = 0 ]; then
        echo "-------------------- COARRAY SMOKE TEST --------------------"
        ctest -R '^coarrays$' --no-tests=error --output-on-failure || GATE_RC=$?
    fi
fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
