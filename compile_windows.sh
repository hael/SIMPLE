#!/bin/bash
# Tests are built by default (simple_test_exec and the *_tester modules) and the fast
# gate runs before installation; --exclude-tests builds the library and executables only.
BUILD_TESTS=ON
ROOT="$(cd "$(dirname "$0")" && pwd)"
source "$ROOT/scripts/simple_build_jobs.sh" || exit $?
for arg in "$@"; do
    case "$arg" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        -h|--help)
            echo "usage: $(basename "$0") [--exclude-tests]"
            echo "Default jobs: available CPUs, or CMAKE_BUILD_PARALLEL_LEVEL when set."
            exit 0
            ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
JOBS=$(simple_build_jobs) || exit $?
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} -G "Unix Makefiles" ..
make -j"$JOBS" || exit $?
# Unless --exclude-tests is given, the build-time test gate (scripts/run_fast_gate.sh)
# runs between build and install: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make -j"$JOBS" install || exit $?; }
exit ${GATE_RC:-0}
