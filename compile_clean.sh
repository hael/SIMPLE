#!/bin/bash
# Clean release build and install into build/.
#
#   ./compile_clean.sh                  library and executables only (default)
#   ./compile_clean.sh --compile-tests  also build simple_test_exec and the test
#                                       programs in production/tests
#
# The test code adds ~210 s of compile CPU and one static link per program; it
# is rarely needed for everyday work, so it is off unless asked for.
BUILD_TESTS=OFF
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --compile-tests) BUILD_TESTS=ON ;;
        -h|--help) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "compile_clean.sh: unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake .. -DBUILD_TESTS=${BUILD_TESTS}
make -j || exit $?
# With --compile-tests, the build-time test gate (scripts/run_fast_gate.sh) runs
# between build and install, as in X: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
