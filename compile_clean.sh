#!/bin/bash
# Clean release build and install into build/.
#
#   ./compile_clean.sh                  library, executables and tests; the fast
#                                       test gate runs before installation (default)
#   ./compile_clean.sh --exclude-tests  library and executables only, no gate
# High-level tests are registered but never run here; invoke them explicitly
# after building with: (cd build && ctest -L highlevel --output-on-failure)
#
# The test code adds ~210 s of compile CPU; --exclude-tests skips it when only
# the executables are needed.
BUILD_TESTS=ON
ROOT="$(cd "$(dirname "$0")" && pwd)"
for arg in "$@"; do
    case "$arg" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        -h|--help) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "compile_clean.sh: unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
rm -rf build
mkdir build
cd build
cmake .. -DBUILD_TESTS=${BUILD_TESTS}
make -j || exit $?
# Unless --exclude-tests is given, the build-time test gate (scripts/run_fast_gate.sh)
# runs between build and install, as in X: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make install || exit $?; }
exit ${GATE_RC:-0}
